# Packing Circles in a Circle -- Solving for minimum radius
from pyomo.environ import *
import numpy as np
import time
import matplotlib.pyplot as plt


#####################################################################
#######################  F U N C T I O N S ##########################
#####################################################################


def plotCylinder():
	return


def plotSphere(fig,ax,x_c,y_c,z_c):
	u, v = np.mgrid[0:2 * np.pi:20j, 0:np.pi:20j]
	x = x_c + np.cos(u) * np.sin(v)
	y = y_c + np.sin(u) * np.sin(v)
	z = z_c + np.cos(v)
	#ax.plot_surface(x, y, z, cmap=plt.cm.YlGnBu_r)
	ax.plot_surface(x, y, z)
	return fig, ax


def plotConfiguration(model):


	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	#ax.set_aspect('equal')
	for ii in model.I:
		# Grab the solutions
		x_ii = value(model.x[ii])
		y_ii = value(model.y[ii])
		z_ii = value(model.z[ii])

		fig,ax = plotSphere(fig,ax,x_c=x_ii,y_c=y_ii,z_c=z_ii)

	#ax.set_box_aspect(aspect=(1,1,1))
	plt.show()

	return

def packingDistance(model,ii,jj):
	x_ii, y_ii, z_ii, w_ii = model.x[ii], model.y[ii], model.z[ii], model.w[ii]
	x_jj, y_jj, z_jj, w_jj = model.x[jj], model.y[jj], model.z[jj], model.w[jj]
	
	distance = ( (x_ii**2 + y_ii**2 + z_ii**2) 
		- 2*w_ii*w_jj*(x_ii*x_jj + y_ii*y_jj + z_ii*z_jj)
		+ (x_jj**2 + y_jj**2 + z_jj**2) ) >= w_ii*w_jj*4	
	return distance


def spherePack(n,radius,height):
	''' 
	Writing a MINLP for packing spheres of unit radius in a cylinder
	'''

	model = ConcreteModel()

	model.N = n 					# An integer
	model.I = RangeSet(1,model.N)	# Indices for Variables

	model.x = Var(model.I,domain=Reals,bounds=(-radius,radius))
	model.y = Var(model.I,domain=Reals,bounds=(-radius,radius))
	model.z = Var(model.I,domain=NonNegativeReals,bounds=(0,height))

	model.w = Var(model.I,domain=Binary)

	model.boundaries = ConstraintList()
	model.overlap = ConstraintList()

	for ii in model.I:
		x_ii, y_ii, z_ii, w_ii = model.x[ii], model.y[ii], model.z[ii], model.w[ii]

		
		# Circle Center Constraint
		model.boundaries.add(expr= (x_ii**2 + y_ii**2) <= (radius - 1)**2)

		model.boundaries.add(expr= z_ii <= w_ii*(height - 1))
		model.boundaries.add(expr= z_ii >= w_ii* 1 )

		model.boundaries.add(expr= x_ii <= w_ii*(radius-1))
		model.boundaries.add(expr= x_ii >= -w_ii*(radius-1))

		model.boundaries.add(expr= y_ii <= w_ii*(radius-1))
		model.boundaries.add(expr= y_ii >= -w_ii*(radius-1))

		# Adding Overlap Constraints
		temp_range = RangeSet(ii+1,model.N)
		for jj in temp_range:
			model.overlap.add(expr= packingDistance(model,ii,jj))

	# quicksum builds the model quicker than sum()
	model.obj = Objective(expr=quicksum(model.w[ii] for ii in model.I),sense=maximize) 

	return model




def packingAlgo(cylinder_diameter,cylinder_height,max_spheres):


	cylinder_radius = 0.5*cylinder_diameter

	# initiating theses variables
	current_best = 0
	last_best = 0
	current_elapsed	= 5 # seconds

	solver = SolverFactory('couenne')

	for ii in range(1,max_spheres+1):

		last_best = current_best

		temp_model = spherePack(n=ii,radius=cylinder_radius,height=cylinder_height)

		# fixing binary variables to reduce load on solver
		for jj in range(1,last_best+1):

			temp_model.w[jj].fix(1)

			# Warm starting model
			temp_model.x[jj] = value(last_model.x[jj])
			temp_model.y[jj] = value(last_model.y[jj])
			temp_model.z[jj] = value(last_model.z[jj])

		tic = time.time()
		try:
			# Try to solve the problem within the alloted time limit
			solver.solve(temp_model, timelimit= 10*current_elapsed) # Hard coded here
			toc = time.time()
			current_elapsed = toc-tic	# seconds

			current_best = int( value(temp_model.obj) )		# storing current best soln
			last_model = temp_model		# recording last soln

			print('----------')
			print('Packing ', ii, ' Spheres')
			print('Time Elapsed is %2.2f seconds' %current_elapsed)


			# If the solver above cannot pack anymore within the time limit
			if current_best == last_best:
				break

		except Exception:
			# The above will try to fit as many as possible and it will timeout for
			# the last one
			toc = time.time()
			current_elapsed = toc-tic	# seconds

			break

	print('----------')
	print('Failed Packing of ', current_best+1, 'spheres')
	print('Time Elapsed is %2.2f seconds' %current_elapsed)
	print('----------')
	print('Maximum Number of Packed Spheres is ', last_best)

	packing_n = last_best
	best_model = last_model
	return best_model, packing_n


#####################################################################
#######################  M A I N - P R O G R A M ####################
#####################################################################

# Cylinder Properties
cylinder_diameter = 10 	
cylinder_radius = 0.5 * cylinder_diameter
cylinder_height = 2

unit_diameter = 2 # fixed constant -- DO NOT CHANGE

cylinder_volume = np.pi * cylinder_radius**2 * cylinder_height
unit_sphere_volume = np.pi *(4/3)

print('Cylinder Volume is ', cylinder_volume)
print('Marble Volume is ', unit_sphere_volume)

estim_packing_n = cylinder_volume / unit_sphere_volume
max_spheres = int(np.ceil(estim_packing_n))

total_tic = time.time()
packing_model, actual_packing_n = packingAlgo(cylinder_diameter,cylinder_height,max_spheres)
total_toc = time.time()

total_elapsed = total_toc - total_tic
packing_factor = actual_packing_n / estim_packing_n

print('----------')
print('Total Time Elapsed is %2.2f seconds' %total_elapsed)
print('Packing Factor is %2.2f' %packing_factor)

plotConfiguration(model=packing_model)