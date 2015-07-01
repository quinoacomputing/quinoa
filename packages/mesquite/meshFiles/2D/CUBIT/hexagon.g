! Database Title                             exo2txt
cubit(hexagon.g): 12/14/2001: 16:04:17
! Database initial variables
         3      2.05               ! dimensions, version number
         7        6         1     ! nodes, elements, element blocks
         0         0               ! #node sets, #side sets
         0         0               ! len: node set list, dist fact length
         0         0         0     ! side sets len: element, node , dist
fact
! Coordinate names
x                                y                                z
! Coordinates
 0.4   0.3 0.000
 1.000 0.000 0.000
 0.500 0.8660254 0.00
 -.500 0.8660254 0.00
 -1.00 0.000 0.00
 -.500 -.8660254 0.00
 0.500 -.8660254 0.00
 
! Node number map
sequence 1..numnp
! Element number map
sequence 1..numel
! Element order map
sequence 1..numel
! Element block    1
         1        6      TRI3      ! ID, elements, name
         3         0      ! nodes per element, attributes
! Connectivity
      1 2 3
      1 3 4
      1 4 5
      1 5 6
      1 6 7
      1 7 2
! Properties
         1            ! Number of ELEMENT BLOCK Properties
! Property Name:
ID
! Property Value(s):
         1
         0            ! Number of NODE SET Properties
         0            ! Number of SIDE SET Properties
! QA Records
         2      ! QA records
CUBIT
7.0b
12/14/2001
16:04:17
exo2txt
 1.13
20011214
16:04:13
! Information Records
         0      ! information records
! Variable names
         0         0         0      ! global, nodal, element variables
