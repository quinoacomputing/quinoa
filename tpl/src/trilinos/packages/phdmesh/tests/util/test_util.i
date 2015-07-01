containers
#
tpi 1
tpi 2
tpi 3
tpi 4
tpi 5
tpi 6
#
# tpi_chunk { target = <Mflops> ; len_array = <#> ; chunk = <#> <#> <#> ... ; }
threadpool 1
tpi_chunk { target = 1000 ; }
#
threadpool 2
tpi_chunk { target = 1000 ; }
#
threadpool 3
tpi_chunk { target = 1000 ; }
#
threadpool 4
tpi_chunk { target = 1000 ; }
#

# comm_all { length =    10000 ; neighbors = -1 1 ; dense = 0 ; }
# comm_all { length =   100000 ; neighbors = -1 1 ; dense = 0 ; }
# comm_all { length =  1000000 ; neighbors = -1 1 ; dense = 0 ; }
# comm_all { length = 10000000 ; neighbors = -1 1 ; dense = 0 ; }
# comm_all { length = 10000000 ; dense = 0 ; }
#
#
bounds
#
sparse
#
dense
#
# global_box <#boxes>
global_box 100
#
oct_tree
#
# oct_tree_part_course <#keys> <#cuts>
oct_tree_part_course 10000 16
#
# oct_tree_comm_part <#keys>
oct_tree_comm_part   10000
#
# oct_tree_global_search <#keys>
oct_tree_global_search 10000
#
oct_tree_global_search_time


