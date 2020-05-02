
/* Maximal edge matching using MPI and UPCXX (RPC) */

#pragma once
#ifndef MAXEMATCHUPXRPC_HPP
#define MAXEMATCHUPXRPC_HPP

#include "graph.hpp"

#include <upcxx/upcxx.hpp>

#include <numeric>
#include <utility>
#include <cstring>
#include <cassert>

// TODO FIXME change comm_world to comm_

class MaxEdgeMatchUPXRPC
{
    public:
        MaxEdgeMatchUPXRPC(Graph* g): 
            g_(g), D_(0), M_(0), mate_(nullptr),
            ghost_count_(nullptr), dobj(this)
        {
            MPI_Comm_size(MPI_COMM_WORLD, &size_);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank_);

            const GraphElem lnv = g_->get_lnv();

            // allocate
            mate_ = new GraphElem[lnv];
            ghost_count_ = new GraphElem[lnv];
            
            // initialize
            std::fill(mate_, mate_ + lnv, -1);
            std::fill(ghost_count_, ghost_count_ + lnv, 0);

            // populate counter that tracks number
            // of ghosts not owned by me
            // ghosts owned per process

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
                        ghost_count_[i] += 1;
                    }
                }
            }
        }

        ~MaxEdgeMatchUPXRPC() {}

        void clear()
        {
            M_.clear();
            D_.clear();
                     
            delete []ghost_count_;
            delete []mate_;
        }
       
        // TODO FIXME not expecting a, b to
        // be large, if large then following
        // absolute tolerance test will fail:
        // http://realtimecollisiondetection.net/blog/?p=89
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
                M_global = new GraphElem[m_global_size];

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

            GraphElem* M_global = nullptr;
            int* rcounts = nullptr;
            int* rdispls = nullptr;

            // communication params
            if (rank_ == 0)
            {
                rcounts = new int[size_];
                rdispls = new int[size_];
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
            maxematch_upxrpc();
            return M_;
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
        
        // initiate RPCs
        
        // reject
        void Push_REJECT(int target, GraphElem data[2])
        {
            assert(target < size_);
            GraphElem data_0 = data[0];
            GraphElem data_1 = data[1];
            current = upcxx::when_all(current,
                    upcxx::rpc(target, [data_0, data_1](upcxx::dist_object<MaxEdgeMatchUPXRPC*>& dobj) {
                        MaxEdgeMatchUPXRPC *here = *dobj;

                        here->deactivate_edge(data_0, data_1);

                        // recalculate mate[x]
                        if (here->mate_[here->g_->global_to_local(data_0)] == data_1) {
                            here->find_mate(data_0);
                        }
                        }, dobj));
        }
       
        // request
        void Push_REQUEST(int target, GraphElem data[2])
        {
            assert(target < size_);
            GraphElem data_0 = data[0];
            GraphElem data_1 = data[1];

            current = upcxx::when_all(current,
                    upcxx::rpc(target,
                        [data_0, data_1](upcxx::dist_object<MaxEdgeMatchUPXRPC*>& dobj) {
                        bool matched = false;
                        MaxEdgeMatchUPXRPC *here = *dobj;

                        // check if y is already matched
                        if (!here->is_matched(data_0))
                        { 
                            // deactivate edge
                            here->deactivate_edge(data_0, data_1);

                            if (here->mate_[here->g_->global_to_local(data_0)] == data_1)
                            {
                                here->M_.emplace_back(data_0, data_1, 0.0);
                                here->D_.push_back(data_0);
                                here->D_.push_back(data_1);

                                matched = true;
                            }
                        } 
                        // put REJECT if matching not possible
                        if (!matched)
                        {
                            // deactivate edge
                            here->deactivate_edge(data_0, data_1);

                            GraphElem sdata[2] = {data_1, data_0}; // {g,l}
                            const int source = here->g_->get_owner(sdata[0]);
                            here->Push_REJECT(source, sdata);
                        }
                    }, dobj));
        }

        // invalid
        void Push_INVALID(int target, GraphElem data[2])
        {
            GraphElem data_0 = data[0];
            GraphElem data_1 = data[1];
            current = upcxx::when_all(current, upcxx::rpc(target,
                    [data_0, data_1](upcxx::dist_object<MaxEdgeMatchUPXRPC*>& dobj) {
                        MaxEdgeMatchUPXRPC *here = *dobj;
                        here->deactivate_edge(data_0, data_1);
                    }, dobj));
        }

        // maximal edge matching using MPI RMA
        void maxematch_upxrpc()
        {
            /* Part #1 -- initialize candidate mate */
            current = upcxx::make_future();

            GraphElem lnv = g_->get_lnv();
            for (GraphElem i = 0; i < lnv; i++)
                find_mate(g_->local_to_global(i));

            /* Part 2 -- compute mate */

            while(1)
            {
                prev = current;
                current = upcxx::make_future();
                prev.wait();

                upcxx::barrier();
                
                do_matching();

                // exit criteria: check if all cross edges have been processed
                GraphElem count = std::accumulate(ghost_count_, ghost_count_ + lnv, 0);
                //std::cout << "count: " << count << std::endl;     
                
                GraphElem total_count = upcxx::reduce_all(count,
                        std::plus<GraphElem>()).wait();

                if (total_count == 0)
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
                    Push_REQUEST(y_owner, y_x);
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
                            Push_INVALID(z_owner, z_x);  
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
                            Push_REJECT(x_owner, x_v);
                        }
                    }
                }   
            }
        }

    private:
        Graph* g_;
        /*
         * List of vertices we've found a match for and whose neighborhood we
         * now want to process. Is inserted in to during each iteration and then
         * emptied.
         */
        std::vector<GraphElem> D_;

        // Stores the list of edges in the match
        std::vector<EdgeTuple> M_;

        GraphElem *mate_;
        // default constructed futures that will
        // be move assigned
        upcxx::future<> current;
        upcxx::future<> prev;

        // count of ghost vertices not owned by me
        GraphElem *ghost_count_; 
        std::vector<int> targets_; // neighbor processes

        int rank_, size_;

        upcxx::dist_object<MaxEdgeMatchUPXRPC*> dobj;
};

#endif
