FLAGS=-O3
# -xmic-avx512 -DCHECK_RESULTS
all: match_rpc match_rma

match_rpc: main.cpp maxematch_upx_rpc.hpp
	upcxx $(FLAGS) -DUSE_UPX_RPC main.cpp -o $@

match_rma: main.cpp maxematch_upx.hpp
	upcxx $(FLAGS) main.cpp -o $@

clean:
	rm -rf match_rpc match_rma
