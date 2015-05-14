import math
import numpy as np
from numpy.linalg import inv
import time
from time import sleep
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#f = open('pos.txt','w') #open file

def inertiaMatrix(pos):
	pos1 = pos.item(0,0)
	pos2 = pos.item(1,0)
	pos3 = pos.item(2,0)
	m11 = 1195000*math.cos(2*pos2+pos3)+970500*math.cos(2*pos2)+639200*math.cos(2*pos2+2*pos3)+620400*math.cos(pos2+pos3)+906300*math.cos(pos2)+1195000*math.cos(pos3)+2189000
	m12 = 3567
	m13 = 3567
	m21 = 3567
	m22 = 4780000*math.cos(pos3/2)*math.cos(pos3/2)+3155000
	m23 = 2390000*math.cos(pos3/2)*math.cos(pos3/2)+1248000
	m31 = 3567
	m32 = 2390000*math.cos(pos3/2)*math.cos(pos3/2)+1248000
	m33 = 2443000

	IM = np.matrix([[m11,m12,m13],[m21,m22,m23],[m31,m32,m33]])
	return IM

def coriolisMatrix (pos,vel):
	pos1 = pos.item(0,0)
	pos2 = pos.item(1,0)
	pos3 = pos.item(2,0)
	vel1 = vel.item(0,0)
	vel2 = vel.item(1,0)
	vel3 = vel.item(2,0)
	c11 = -vel2*(1195000*math.sin(2*pos2+pos3)+970500*math.sin(2*pos2)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+453100*math.sin(pos2))-vel3*(597500*math.sin(2*pos2+pos3)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+597500*math.sin(pos3))
	c12 = -vel1*(1195000*math.sin(2*pos2+pos3)+970500*math.sin(2*pos2)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+453100*math.sin(pos2))
	c13 = -vel1*(597500*math.sin(2*pos2+pos3)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+597500*math.sin(pos3))
	c21 = vel1*(1195000*math.sin(2*pos2+pos3)+970500*math.sin(2*pos2)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+453100*math.sin(pos2))-vel2*(1195000*math.sin(2*pos2+pos3)+970500*math.sin(2*pos2)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+453100*math.sin(pos2))-vel2*(((4643513087359461*math.cos(pos2)*math.sin(pos2))/618970019642690137449562112)+(1242272829422879*math.cos(pos2+pos3)*math.sin(pos2+pos3))/154742504910672534362390528)-vel3*((597500*math.sin(2*pos2+pos3)+639200*math.sin(2*pos2+2*pos3)+310200*math.sin(pos2+pos3)+597500*math.sin(pos3))+(1242272829422879*math.cos(pos2+pos3)*math.sin(pos2+pos3))/154742504910672534362390528)
	c22 = 2390000*vel3*math.cos(pos3/2)*math.sin(pos3/2)
	c23 = -2390000*vel2*math.cos(pos3/2)*math.sin(pos3/2) - 2390000*vel3*math.cos(pos3/2)*math.sin(pos3/2)
	c31 = -vel2*(5481523169798827*math.cos(pos2+pos3)*math.sin(pos2+pos3)/11150372599265311570767859136324180752990208 + 2390000*math.cos(pos3/2)*math.sin(pos3/2)) - vel3*((5481523169798827*math.cos(pos2+pos3)*math.sin(pos2+pos3)/11150372599265311570767859136324180752990208) + 2390000*math.cos(pos3/2)*math.sin(pos3/2))
	c32 = 2390000*vel2*math.cos(pos3/2)*math.sin(pos3/2)
	c33 = 0
	CM = np.matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])
	return CM

def gravityMatrix(pos):
	pos1 = pos.item(0,0)
	pos2 = pos.item(1,0)
	pos3 = pos.item(2,0)
	g1 = 0
	g2 = 937760463*math.cos(pos2+pos3)/10000 + 412323129*math.cos(pos2)/2500
	g3 = 937760463*math.cos(pos2+pos3)/10000
	GM = np.matrix([[g1],[g2],[g3]])
	return GM

def FK (pos): #Formard Kinematics
	a1 = 22
	a2 = 125
	a3 = 150
	alpha1=math.pi/2
	alpha2=0
	alpha3=0
	d1 = 71
	d2 = 0
	d3 = 0
	theta1 = pos.item(0,0)
	theta2 = pos.item(1,0)
	theta3 = pos.item(2,0)
	A1 = [[np.cos(theta1), -np.cos(alpha1)*np.sin(theta1), np.sin(alpha1)*np.sin(theta1), a1*np.cos(theta1)], [np.sin(theta1),np.cos(alpha1)*np.cos(theta1),-np.sin(alpha1)*np.cos(theta1),a1*np.sin(theta1)],[0, np.sin(alpha1),np.cos(alpha1),d1],[0,0,0,1]]
	A2 = [[np.cos(theta2), -np.cos(alpha2)*np.sin(theta2), np.sin(alpha2)*np.sin(theta2), a2*np.cos(theta2)], [np.sin(theta2),np.cos(alpha2)*np.cos(theta2),-np.sin(alpha2)*np.cos(theta2),a2*np.sin(theta2)],[0, np.sin(alpha2),np.cos(alpha2),d2],[0,0,0,1]]
	A3 = [[np.cos(theta3), -np.cos(alpha3)*np.sin(theta3), np.sin(alpha3)*np.sin(theta3), a3*np.cos(theta3)], [np.sin(theta3),np.cos(alpha3)*np.cos(theta3),-np.sin(alpha3)*np.cos(theta3),a3*np.sin(theta3)],[0, np.sin(alpha3),np.cos(alpha3),d3],[0,0,0,1]]
	T03_int = np.dot(A1,A2)
	T03 = np.dot(T03_int,A3)

	x = T03.item((0,3))+30*np.cos(theta1)
	y = T03.item((1,3))+30*np.sin(theta1)
	z = T03.item((2,3))+70

	R = np.matrix([[T03.item(0,0),T03.item(0,1),T03.item(0,2)],[T03.item(1,0),T03.item(1,1),T03.item(1,2)],[T03.item(2,0),T03.item(2,1),T03.item(2,2)]])
	F = np.matrix([[x],[y],[z]])
	
	simpleX = np.cos(theta1)*(52+150*np.cos(theta3-theta2)+125*np.cos(theta2));
	simpleY = np.sin(theta1)*(52+150*np.cos(theta3-theta2)+125*np.cos(theta2));
	simpleZ = 125*np.sin(theta2)+150*np.sin(theta2-theta3)+141;
	simple = np.matrix([[simpleX],[simpleY],[simpleZ]])
	return simple

#General Parameters for Faulhaber 1524006SR with a 141:1 gear head
Ke = 0.628	#Back EMF Constant 
Km = 6.0 	#Torque Constant
R = 5.1  	#Terminal Resistance
B = Ke*Km/R
xi = 1.0  	#Damping factor (for critical damping system, xi = 1)
w = 10.0 	#Hz
J = 0.58	#Motor inertia

#Kp and Kd
Kp = w*w*J
Kd = 2*xi*w*J-B
Ki = 0

#initialize data matrices
pos = np.matrix([[0.0],[0.0],[0.0]])
vel = np.matrix([[0.0],[0.0],[0.0]])
acc = np.matrix([[0.0],[0.0],[0.0]])
T =  np.matrix([[0.0],[0.0],[0.0]])
accDes = np.matrix([[0.0],[0.0],[0.0]])
velDes = np.matrix([[0.0],[0.0],[0.0]])
posDes = np.matrix([[0.0],[0.0],[0.0]])
posDesPrev = np.matrix([[0.0],[0.0],[0.0]])
accDesPrev = np.matrix([[0.0],[0.0],[0.0]])
controlInput = np.matrix([[0.0],[0.0],[0.0]])

DT = 0.01 	#initial delta time
i = 0 		#sample number
t = 0 		#Total time

#initialize array for plotting
Array1 = []
Array2 = []
Array3 = []
Array4 = []
Array5 = []
Array6 = []
timeArray = []

#Needed for realtim plotting
#plt.ion()

while True:
	posDes = np.matrix([[math.sin(t)],[math.cos(t)],[t]]) #Trajectory
	time_initial = time.time()
	IM = inertiaMatrix(pos)
	CM = coriolisMatrix(pos,vel)
	GM = gravityMatrix(pos)
	#T = np.dot(IM,acc) + np.dot(CM,vel) + GM
	acc = np.dot(inv(IM),(T-GM-np.dot(CM,vel)))
	
	#Integral of Acceleration -> velocity
	acc_int = acc*DT
	vel = acc_int + vel

	#Integral of Velocity -> position
	vel_int = vel*DT
	pos = vel_int + pos

	#Derivative of desired position, posDes -> velDes
	velDes = (posDes - posDesPrev)/DT
	posDesPrev = posDes

	#Derivative of desired velocity, velDes -> accDes
	accDes = (accDes - accDesPrev)/DT
	accDesPrev = accDes

	#Velocity error (ed')
	velErr = velDes - vel

	#Desired position error (ed'')
	posErr = posDes - pos

	#Update forward kinematics
	fk = FK(pos)

	#Plotting arrays for forward kinematics
	Array1.append(pos.item(0,0)) #-> X
	Array2.append(pos.item(1,0)) #-> Y
	Array3.append(pos.item(2,0)) #-> Z
	Array4.append(posDes.item(0.0))
	Array5.append(posDes.item(1.0))
	Array6.append(posDes.item(2.0))

	#plt.scatter(Array2,Array1)
	#plt.draw()

	#PD Control
	T = np.dot(IM,(accDes+Kd*velErr+Kp*posErr)) + np.dot(CM,vel) + GM

	print "PosEr1:", pos.item(0,0),"PosEr2:",pos.item(1,0), "PosEr3:", pos.item(2,0)
	#print DT
	#Write to file
	#f.write(str(math.degrees(pos.item(0,0))) + " " + str(math.degrees(pos.item(1,0))) + " " + str(math.degrees(pos.item(2,0))) + "\n")

	#Emergency for DT
	#If DT is small, integral gives error
	
	while (time.time() - time_initial) < 0.001:
		time.sleep(0.001)
	

	DT = time.time() - time_initial	#calculate DT
	t = t + DT 	#update time
	timeArray.append(t)	#update time array for plotting

	if i > 10000:
		break
	i=i+1
print "Total Time:", t


#3D Plotting
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(Array1,Array2,Array3, c='r') #pos
ax.plot(Array4,Array5,Array6, c='b') #pos error
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
#plt.plot(timeArray,Array2)
#plt.plot(timeArray,Array1)
#plt.ylabel('Pos1 (Degrees)')
#plt.xlabel('Time')
plt.show()
#f.close()	#Close file