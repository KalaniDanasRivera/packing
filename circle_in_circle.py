# Packing Circles in a Circle -- Solving for minimum radius
from pyomo.environ import *
import numpy as np
import time
import matplotlib.pyplot as plt



def plotConfiguration(model):
	theta = np.linspace(0,360,100, dtype='float')
	unit_circle = np.array([np.cos(np.deg2rad(theta)).T, np.sin(np.deg2rad(theta)).T])
	outer_circle = value(model.radius()) * unit_circle
	
	
	plt.figure()
	plt.fill(outer_circle[0,:], outer_circle[1,:])


	for ii in model.I:
		temp_center = np.dot(np.diag([value(model.x[ii]) , value(model.y[ii])]), np.ones((2,100)) )
		temp_circle = temp_center + unit_circle

		plt.fill(temp_circle[0,:], temp_circle[1,:],alpha=0.5)

	plt.axis('equal')
	plt.show()
	return


def circlePack(n):
	''' 
	Writing a NLP for packing circles of unit radius in an outer circle
	'''

	bound = 5


	model = ConcreteModel()

	model.N = n 					# An integer
	model.I = RangeSet(1,model.N)	# Indices for Variables

	model.radius = Var(domain=NonNegativeReals,bounds=(1,bound))
	model.x = Var(model.I,domain=Reals,bounds=(-bound,bound))
	model.y = Var(model.I,domain=Reals,bounds=(-bound,bound))

	model.boundaries = ConstraintList()
	model.overlap = ConstraintList()


	model.boundaries.add(expr= model.radius >= model.N**(0.5))

	for ii in model.I:
		x_ii, y_ii = model.x[ii], model.y[ii]

		# Circle Center Constraint
		model.boundaries.add(expr= x_ii**2 + y_ii**2 <= (model.radius - 1)**2)

		temp_range = RangeSet(ii,model.N)

		for jj in temp_range:
			if jj != ii:
				x_jj, y_jj = model.x[jj], model.y[jj]

				delta_x = x_ii - x_jj
				delta_y = y_ii - y_jj

				model.overlap.add(expr= delta_x**2 + delta_y**2 >= 4)
				#model.overlap.add(expr= delta_x**2 + delta_y**2 <= (2*model.radius-2)**2)

	model.obj = Objective(expr=model.radius)
	return model

max_circles = 5

for ii in range(1,max_circles+1):

	temp_model = circlePack(n=ii)

	tic = time.time()
	SolverFactory('couenne').solve(temp_model)
	toc = time.time()

	print('----------')
	print('Number of circles: ', ii)
	print('Radius is ', value(temp_model.radius))
	print('Time Elapsed is ', toc-tic, ' seconds')
	print('----------')


#plotConfiguration(temp_model)
