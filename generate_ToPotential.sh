#!/bin/bash

out_file_name='ToPotential'

grid_step=0.05
grid_max=5.0

rm -f ${out_file_name}

# to simulate imaginary-theta term with multicanonic, choose topo-potetial V(i) = -im_theta*i
im_theta=0.5

# print topo-potential to file
for i in $( seq -${grid_max} ${grid_step} ${grid_max} ); do
	# substitute here the desired topo-potential as a function of the grid dummy variable i
	#							|______________________________________________________
	#                                                                    |
	#                                                                    v
	awk -v i=${i} -v t=${im_theta} 'BEGIN{ printf "%.5lf %.18lf\n", i, -t*i  }' >> ${out_file_name}
done
