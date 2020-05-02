
/* Maximal edge matching using MPI and UPCXX */

#pragma once
#ifndef MAXEMATCHUPX_HPP
#define MAXEMATCHUPX_HPP

#include "graph.hpp"

#include <upcxx/upcxx.hpp>

#include <numeric>
#include <utility>
#include <cstring>
#include <cassert>

#define MATE_REQUEST        1 // mate[x] = y, if mate[y] = x, then match (x/y in different processes)
#define MATE_REJECT         2 // reject vertex mate, requires to update mate
#define MATE_INVALID        3 // invalidate edge

// TODO FIXME change comm_world to comm_

class MaxEdgeMatchUPX
{
    public:
        MaxEdgeMatchUPX(Graph* g): 
            g_(g), D_(0), M_(0), mate_(nullptr),
            qbuf_(nullptr), ghost_count_(nullptr),
            nghosts_(0), nghosts_indices_(0), rdispls_(0),
            scounts_(NULL), rcounts_(NULL), prcounts_(NULL),
            wbuf_(nullptr), gptrs_(0), futs_(upcxx::make_future())
        {
            MPI_Comm_size(MPI_COMM_WORLD, &size_);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

            const GraphElem lnv = g_->get_lnv();

            // allocate
            mate_ = new GraphElem[lnv];
            assert(mate_);
            ghost_count_ = new GraphElem[lnv];
            assert(ghost_count_);
            
            // initialize
            std::fill(mate_, mate_ + lnv, -1);
            std::fill(ghost_count_, ghost_count_ + lnv, 0);

            rdispls_.resize(size_, 0); // stores remote displacement
            nghosts_.resize(size_, 0);
            nghosts_indices_.resize(size_, 0);
            
            // populate counter that tracks number
            // of ghosts not owned by me
            // ghosts owned per process
            GraphElem tot_ghosts = 0;

            for (GraphElem i = 0; i < lnv; i++)
            {
                GraphElem e0, e1;
                g_->edge_range(i, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    Edge const& edge = g_->get_edge(e);

                    // all edge is active in the beginning,
                    // so no need to check edge.active_
                    const int p = g_->get_owner(edge.tail_);
                    if (p != rank_)
                    {
                        if (std::find(targets_.begin(), targets_.end(), p) 
                                == targets_.end())
                            targets_.push_back(p);
                        nghosts_[p] += 1;
                        ghost_count_[i] += 1;
                    }
                }

                tot_ghosts += ghost_count_[i];
            }

            // allocate outgoing buffer
            // 2 because we can process a vertex at most twice
            // 3 because we send 3 data parameters
            nelems_ = tot_ghosts*2*3;
            qbuf_ = new GraphElem[nelems_];
            assert(qbuf_);


            // prefix sum for calculating 
            // indices of outgoing buffers
             
            // scan to compute remote 
            // displacements for RMA CALLS
            GraphElem disp = 0;
            for (int p = 0; p < size_; p++)
            {
                nghosts_indices_[p] = disp;
                disp += nghosts_[p]*2*3;
            }
               
            // incoming updated prefix sums
            MPI_Alltoall(nghosts_indices_.data(), 1, MPI_GRAPH_TYPE, 
                    rdispls_.data(), 1, MPI_GRAPH_TYPE, MPI_COMM_WORLD);
           
            scounts_ = new GraphElem[size_];
            assert(scounts_);
            rcounts_ = new GraphElem[size_];
            assert(rcounts_);
            prcounts_ = new GraphElem[size_];
            assert(prcounts_);
            for (int i = 0; i < size_; i++) {
                scounts_[i] = 0;
                rcounts_[i] = 0;
                prcounts_[i] = 0;
            }
        }
        
        // every process creates global_ptr
        void create_upx_ptr(upcxx::dist_object<upcxx::global_ptr<GraphElem>>& dobj)
        {
            for (int p = 0; p < targets_.size(); p++)
            {
                gptrs_.insert(std::make_pair< GraphElem, upcxx::global_ptr<GraphElem> >
                        (targets_[p], dobj.fetch(targets_[p]).wait())); 
            }
            upcxx::barrier();
            
            wbuf_ = dobj->local();
            assert(dobj->is_local());
        }

        // TODO FIXME is there a memory leak?
        void destroy_upx_ptr() 
        { gptrs_.clear(); }

        GraphElem get_nelems() const { return nelems_; }

        ~MaxEdgeMatchUPX() {}

        void clear()
        {
            M_.clear();
            D_.clear();
            
            delete []qbuf_;
            delete []ghost_count_;
            delete []mate_;

            rdispls_.clear();
            delete scounts_;
            delete rcounts_;
            delete prcounts_;
            nghosts_.clear();
            nghosts_indices_.clear();
        }

        inline bool is_same(GraphWeight a, GraphWeight b) 
        { return std::abs(a - b) <= std::numeric_limits<GraphWeight>::epsilon(); }

        /* Validation */
        // if mate[mate[v]] == v then
        // we're good
        void check_results()
        {
            // gather M_ and mate_
            const int lnv = g_->get_lnv();
            unsigned int m_size = M_.size(), m_global_size = 0;
            // i,j
            m_size *= 2;
            GraphElem* M_buf = new GraphElem[m_size];
            assert(M_buf);

            GraphElem* M_global = nullptr;
            GraphElem* mate_global = nullptr;
            
            // communication params from M_ and mate_
            int* rcounts = nullptr;
            int* rdispls = nullptr;
            int* m_rcounts = nullptr;
            int* m_rdispls = nullptr;
            
            // communication params for M
            if (rank_ == 0)
            {
                rcounts = new int[size_];
                rdispls = new int[size_];
                m_rcounts = new int[size_];
                m_rdispls = new int[size_];
                assert(rcounts && rdispls && m_rcounts && m_rdispls);
            }

            // put M_ into a contiguous buffer
            for (int i = 0, j = 0; i < m_size; i+=2, j++)
            {
                M_buf[i]    = M_[j].ij_[0];
                M_buf[i+1]  = M_[j].ij_[1];
            }

            MPI_Gather(&m_size, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Gather(&lnv, 1, MPI_INT, m_rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Reduce(&m_size, &m_global_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            
            // communication params (at root)
            if (rank_ == 0)
            {
                const GraphElem nv = g_->get_nv();
                mate_global = new GraphElem[nv];
                assert(mate_global);
                M_global = new GraphElem[m_global_size];
                assert(M_global);

                unsigned int index = 0, m_index = 0;
                for (int p = 0; p < size_; p++)
                {
                    rdispls[p] = index;
                    index += rcounts[p];
                    m_rdispls[p] = m_index;
                    m_index += m_rcounts[p];
                }
            }
            
            MPI_Barrier(MPI_COMM_WORLD);

            // M_
            MPI_Gatherv(M_buf, m_size, MPI_GRAPH_TYPE, M_global, rcounts, rdispls, 
                    MPI_GRAPH_TYPE, 0, MPI_COMM_WORLD);
            // mate
            MPI_Gatherv(mate_, lnv, MPI_GRAPH_TYPE, mate_global, m_rcounts, m_rdispls, 
                    MPI_GRAPH_TYPE, 0, MPI_COMM_WORLD);
            
            MPI_Barrier(MPI_COMM_WORLD);

            // data gathered, now validate
            if (rank_ == 0)
            {
                bool success = true;
                for (int i = 0; i < m_global_size; i+=2)
                {
                    if ((mate_global[mate_global[M_global[i]]] != M_global[i])
                            || (mate_global[mate_global[M_global[i+1]]] != M_global[i+1]))
                    {
                        std::cout << "Validation FAILED." << std::endl; 
                        std::cout << "mate_[mate_[" << M_global[i] << "]] != " << M_global[i] << " OR " 
                            << "mate_[mate_[" << M_global[i+1] << "]] != " << M_global[i+1] << std::endl;
                        success = false;
                        break;
                    }
                }
                if (success) 
                    std::cout << "Validation SUCCESS." << std::endl;
            }
            
            MPI_Barrier(MPI_COMM_WORLD);

            // clear buffers
            delete []M_global;
            delete []mate_global;
            delete []M_buf;

            delete []rcounts;
            delete []rdispls;
            delete []m_rcounts;
            delete []m_rdispls;
        }

        // print the contents of M_
        void print_M() const
        {
            // gather M_
            unsigned int m_size = M_.size(), m_global_size = 0;
            // i,j
            m_size *= 2;
            GraphElem* M_buf = new GraphElem[m_size];
            assert(M_buf);

            GraphElem* M_global = nullptr;
            int* rcounts = nullptr;
            int* rdispls = nullptr;

            // communication params
            if (rank_ == 0)
            {
                rcounts = new int[size_];
                rdispls = new int[size_];
                assert(rcounts && rdispls);
            }

            // put M_ into a contiguous buffer
            for (int i = 0, j = 0; i < m_size; i+=2, j++)
            {
                M_buf[i]    = M_[j].ij_[0];
                M_buf[i+1]  = M_[j].ij_[1];
            }

            MPI_Gather(&m_size, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Reduce(&m_size, &m_global_size, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

            // communication params (at root)
            if (rank_ == 0)
            {
                M_global = new GraphElem[m_global_size];
                assert(M_global);

                unsigned int index = 0;
                for (int p = 0; p < size_; p++)
                {
                    rdispls[p] = index;
                    index += rcounts[p];
                }
            }

            MPI_Gatherv(M_buf, m_size, MPI_GRAPH_TYPE, M_global, rcounts, rdispls, 
                    MPI_GRAPH_TYPE, 0, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);

            // print mates
            if (rank_ == 0)
            {
                std::cout << "Matched vertices: " << std::endl;
                for (int i = 0; i < m_global_size; i+=2)
                    std::cout << M_global[i] << " ---- " << M_global[i+1] << std::endl;
            }
            
            MPI_Barrier(MPI_COMM_WORLD);

            // clear buffers
            delete []M_global;
            delete []M_buf;
            delete []rcounts;
            delete []rdispls;
        }
                
        /* Maximal edge matching */
        std::vector<EdgeTuple> const& operator()()
        {
            maxematch_upx();
            return M_;
        }

        void NbPut(int tag, int target, GraphElem data[2])
        {
            const GraphElem index = nghosts_indices_[target] + scounts_[target];

            qbuf_[index] = data[0];
            qbuf_[index + 1] = data[1];
            qbuf_[index + 2] = tag;

            // get displacement
            GraphElem tdisp = rdispls_[target] + scounts_[target];

            futs_ = upcxx::when_all(futs_,
                    upcxx::rput(&qbuf_[index], gptrs_[target] + tdisp, 3));
            
            scounts_[target] += 3;
        }
        
        // search v in M_ (local list
        // of matched vertices)
        bool is_matched(const GraphElem v)
        {
            auto found = std::find_if(M_.begin(), M_.end(), 
                    [&](EdgeTuple const& et) 
                    { return ((et.ij_[0] == v) || (et.ij_[1] == v)); });
            if (found == std::end(M_))
                return false;

            return true;
        }

        // process matched vertices
        // in my local queue, ignore
        // ghost vertices
        void do_matching()
        {
            while(!D_.empty())                               
            {    
                GraphElem v = D_.back();
                D_.pop_back();
                const int v_owner = g_->get_owner(v);

                if (v_owner == rank_)
                    process_neighbors(v);
            }
        }
        
        // check matching requests
        void process_window()
        {
            // incoming data sizes
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Alltoall(scounts_, 1, MPI_GRAPH_TYPE, 
                    rcounts_, 1, MPI_GRAPH_TYPE, MPI_COMM_WORLD);
            MPI_Barrier(MPI_COMM_WORLD);
            // access local window and process data
            for (int k = 0; k < size_; k++)
            {
                const GraphElem index = nghosts_indices_[k];
                const GraphElem start = prcounts_[k];
                const GraphElem end = rcounts_[k];
                
                for (int i = start; i < end; i+=3)
                {
                    // REQUEST: may result in a match
                    if (wbuf_[index+i+2] == MATE_REQUEST) 
                    {
                        bool matched = false;

                        // check if y is already matched
                        if (!is_matched(wbuf_[index+i]))
                        { 
                            // deactivate edge
                            deactivate_edge(wbuf_[index+i], wbuf_[index+i+1]);

                            if (mate_[g_->global_to_local(wbuf_[index+i])] 
                                    == wbuf_[index+i+1])
                            {
                                M_.emplace_back(wbuf_[index+i], wbuf_[index+i+1], 0.0);

                                D_.push_back(wbuf_[index+i]);
                                D_.push_back(wbuf_[index+i+1]);

                                matched = true;
                            }
                        } 

                        // put REJECT if matching not possible
                        if (!matched)
                        {
                            // deactivate edge
                            deactivate_edge(wbuf_[index+i], wbuf_[index+i+1]);

                            GraphElem data[2] = {wbuf_[index+i+1], wbuf_[index+i]}; // {g,l}
                            const int source = g_->get_owner(data[0]);

                            NbPut(MATE_REJECT, source, data);
                        }                
                    } 
                    else if (wbuf_[index+i+2] == MATE_REJECT)
                    {
                        deactivate_edge(wbuf_[index+i], wbuf_[index+i+1]);

                        // recalculate mate[x]
                        if (mate_[g_->global_to_local(wbuf_[index+i])] 
                                == wbuf_[index+i+1])
                            find_mate(wbuf_[index+i]);
                    }
                    else if (wbuf_[index+i+2] == MATE_INVALID) { // INVALID: deactivate x -- v
                        deactivate_edge(wbuf_[index+i], wbuf_[index+i+1]);
                    }
                }
                
                // retain past recv counts
                prcounts_[k] = rcounts_[k];
                rcounts_[k] = 0;
            }
        }

        // maximal edge matching using MPI RMA
        void maxematch_upx()
        {
            /* Part #1 -- initialize candidate mate */
            
            GraphElem lnv = g_->get_lnv();
            for (GraphElem i = 0; i < lnv; i++)
                find_mate(g_->local_to_global(i));

            /* Part 2 -- compute mate */
            while(1)
            {
                futs_.wait();
                futs_ = upcxx::make_future();

                process_window();
                do_matching();

                // exit criteria: check if all cross edges have been processed
                GraphElem count = std::accumulate(ghost_count_, ghost_count_ + lnv, 0);
                //std::cout << "count: " << count << std::endl;     

                count = upcxx::reduce_all(count, std::plus<GraphElem>()).wait();
                //MPI_Allreduce(MPI_IN_PLACE, &count, 1, MPI_GRAPH_TYPE, 
                //        MPI_SUM, MPI_COMM_WORLD);

                if (count == 0)
                    break;
            }
        } 

        // expecting v to be local index
        // require global_to_local
        // before passing local computation
        void compute_mate(const GraphElem v, Edge& max_edge)
        {
            GraphElem e0, e1;
            g_->edge_range(v, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.active_)
                {
                    if (edge.edge_.weight_ > max_edge.weight_)
                        max_edge = edge.edge_;

                    // break tie using vertex index
                    if (is_same(edge.edge_.weight_, max_edge.weight_))
                        if (edge.edge_.tail_ > max_edge.tail_)
                            max_edge = edge.edge_;
                }
            }
        }

        // x is owned by me, y may be a ghost
        // deactivate edge x -- y and decrement
        inline void deactivate_edge(GraphElem x, GraphElem y)
        {
            GraphElem e0, e1;
            const GraphElem lx = g_->global_to_local(x);
            const int y_owner = g_->get_owner(y);

            g_->edge_range(lx, e0, e1);

            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                if (edge.edge_.tail_ == y && edge.active_)
                {
                    edge.active_ = false;

                    if (y_owner != rank_)
                        ghost_count_[lx] -= 1;
                    
                    break;
                }
            }
        }

        // x is owned by me
        // compute y = mate[x], if mate[y] = x, match
        // else if y = -1, invalidate all edges adj(x)
        void find_mate(GraphElem x)
        {
            const GraphElem lx = g_->global_to_local(x);
            Edge x_max_edge;

            compute_mate(lx, x_max_edge);
            const GraphElem y = mate_[lx] = x_max_edge.tail_;

            // initiate matching request
            if (y != -1)
            {
                // check if y can be matched
                const int y_owner = g_->get_owner(y);
                if (y_owner == rank_)
                {
                    if (mate_[g_->global_to_local(y)] == x)
                    {
                        D_.push_back(x);
                        D_.push_back(y);
                        M_.emplace_back(x, y, x_max_edge.weight_);

                        // mark y<->x inactive, because its matched
                        deactivate_edge(y, x);
                        deactivate_edge(x, y);
                    }
                }
                else // ghost, initate REQUEST
                {
                    deactivate_edge(x, y);

                    GraphElem y_x[2] = {y, x};
                    NbPut(MATE_REQUEST, y_owner, y_x);  
                }
            }
            else // invalidate all neigboring vertices 
            {
                GraphElem e0, e1;
                g_->edge_range(lx, e0, e1);

                for (GraphElem e = e0; e < e1; e++)
                {
                    EdgeActive& edge = g_->get_active_edge(e);
                    const GraphElem z = edge.edge_.tail_;

                    if (edge.active_)
                    {
                        // invalidate x -- z
                        edge.active_ = false;
                        
                        const int z_owner = g_->get_owner(z);
                        
                        if (z_owner == rank_)
                            deactivate_edge(z, x);
                        else // ghost, initiate INVALID
                        {
                            ghost_count_[lx] -= 1;
                            
                            GraphElem z_x[2] = {z, x};
                            NbPut(MATE_INVALID, z_owner, z_x);  
                        }
                    }
                }
            }
        }

        // process neighbors of matched vertices
        void process_neighbors(GraphElem v)
        {
            GraphElem e0, e1;
            const GraphElem lv = g_->global_to_local(v);
            g_->edge_range(lv, e0, e1);

            // find unmatched vertices in v's neighborhood
            for (GraphElem e = e0; e < e1; e++)
            {
                EdgeActive& edge = g_->get_active_edge(e);
                
                if (edge.active_)
                {
                    const GraphElem x = edge.edge_.tail_;

                    // v is matched with one and only one of its neighbor,
                    // recalculate new mate for all x in N(v) whose candidate 
                    // mate is v
                    // if mate[v] != x then mate[x] != v, for cases where
                    // mate[x] == v, mate[x] needs to be recalculated
                    if (mate_[lv] != x)
                    {
                        // invalidate v - x, because v is already matched, 
                        // and not with x
                        edge.active_ = false;
                        const int x_owner = g_->get_owner(x);

                        // find another mate for x, as v 
                        // is already matched
                        if (x_owner == rank_)
                        {
                            // invalidate x - v
                            deactivate_edge(x, v);

                            // find new candidate
                            if (mate_[g_->global_to_local(x)] == v)
                                find_mate(x);
                        }
                        else // send INVALID to invalidate x-v and recompute mate[x] if necessary
                        {
                            ghost_count_[lv] -= 1; 
                            
                            GraphElem x_v[2] = {x, v};
                            NbPut(MATE_REJECT, x_owner, x_v);
                        }
                    }
                }   
            }
        }

    private:
        Graph* g_;
        std::vector<GraphElem> D_;
        std::vector<EdgeTuple> M_;
        GraphElem *mate_;

        // number of elements on data window
        // and local data buffer
        GraphElem nelems_;
    
        // count of ghost vertices not owned by me
        GraphElem *ghost_count_; 
        
        // outgoing data queue
        GraphElem *qbuf_; 
        GraphElem *countp_;

        // counters
        std::vector<GraphElem> nghosts_, nghosts_indices_, rdispls_;
        GraphElem *scounts_, *rcounts_, *prcounts_;

        // global ptrs to each of my neighbor's data
        std::unordered_map< GraphElem, upcxx::global_ptr<GraphElem> > gptrs_;
        upcxx::future<> futs_;

        std::vector<int> targets_; // neighbor processes
        GraphElem* wbuf_; // ptr to locally owned data

        int rank_, size_;
};

#endif
