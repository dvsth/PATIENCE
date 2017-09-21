import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

#--------------------- fetch simulation settings

# open the config file to read parameters
config = ConfigParser.ConfigParser()
config.read("config.ini")

# field data
B = [ config.getfloat("field" , "B_x") , config.getfloat("field" , "B_y") , config.getfloat("field" , "B_z") ] 

# particle data
mass = config.getfloat("particle" , "mass")
charge = config.getfloat("particle" , "charge")

# boundary conditions
pos_0 = [ config.getfloat("boundary-conditions","x_0") , config.getfloat("boundary-conditions","y_0") , config.getfloat("boundary-conditions","z_0") ]
vel_0 = [ config.getfloat("boundary-conditions","v_x_0") , config.getfloat("boundary-conditions","v_y_0") , config.getfloat("boundary-conditions","v_z_0") ]

#time info
steps = config.getint("time-step","steps")
t_final = config.getfloat("time-step","t_final")
h = t_final/steps

#storage variables for plotting
plot_pos_x = list()
plot_pos_y = list()

#leapfrog variables
vel_plus_half = [ 0. , 0. , 0.]
vel_minus_half = [ 0. , 0. , 0.]

#--------------------- run the simulation

#calculate the magnetic field as a function of time
def field ( t ) :
	return B
	
#calculate the acceleration
def acceleration ( BField , velocity ) :
	return np.multiply(charge/mass , np.cross( velocity , BField ))
	
#calculate the velocity	
def velocityUpdate ( h , v_prev , t ) :
	return np.add(v_prev ,  np.multiply( h , acceleration( field(t) , v_prev )))
	
#calculate the position
def positionUpdate ( h , x_prev , v_prev , v_new , t ) :
	return np.add(x_prev , np.multiply(h/2.0 , np.add(v_new , v_prev)))
	
#algorithm

print ""
print "###########################################"
print "#         commence leapfrogging ...       #"
print "###########################################"
print "#                                         #"               
print "#               _  --  _                  #"
print "#            -            -               #"
print "#         -                  -            #"
print "#       -                      -          #"
print "#      -                        -         #"
print "#     _                          _        #"
print "#                                         #"
print "#          wow, what a leap               #"
print "#                                         #"
print "#  this is one small step for a compiler, #"
print "#     but a giant leap for PATIENCE.      #"
print "#                                         #"  
print "###########################################"
print ""

#compute half-step leapfrog starter
a_0 = np.multiply( charge/mass , np.cross( vel_0 , B ) )
vel_minus_half = np.subtract( vel_0 , np.multiply( h/2. , a_0 ))

#begin integration
for step in range(steps):
	time = step*h
	
	vel_plus_half = velocityUpdate( h , vel_minus_half , time )
	pos_1 = positionUpdate( h , pos_0 , vel_minus_half , vel_plus_half , time )
	
	plot_pos_x.append(pos_1[0])
	plot_pos_y.append(pos_1[1])
	
	vel_minus_half = vel_plus_half
	pos_0 = pos_1
	
#--------------------- generate the plot
fig = plt.figure()

fig.suptitle('Trajectory for a proton \n computed using \n Leapfrog scheme \n in ' + str(steps) + " steps" , fontweight='bold')
pltPos = fig.add_subplot(111)

pltPos.grid(b=True, color='k', linestyle='--')
pltPos.set_xlabel('x-position (meters)' , fontweight='bold')
pltPos.set_ylabel('y-position (meters)' , fontweight='bold')

pltPos.plot(plot_pos_x , plot_pos_y , 'b-', label='path')
pltPos.scatter(plot_pos_x[0] , plot_pos_y[0], color='k' , s=100,marker='o' , label='initial position')
pltPos.scatter(plot_pos_x[-1] , plot_pos_y[-1], color='k' , s=100,marker='x' , label='final position')

pltPos.set_aspect(1.0)
plt.legend(loc='best')
plt.show()



	
	


	
	





