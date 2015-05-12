#PD Position Control for 3-DoF Spatial Robot
#Faulhaber 1524006SR
#Beaglebone Black with Dual Motor Controller Cape (DMCC)

import DMCC # Dual motor controller cape library
import time # Library for calculating the elapsed time between executions
import numpy as np
import math

def trajectory (t):
	x = 10*math.cos(math.pi*t/2*10) + 120
	y = 10*math.sin(math.pi*t/2*10)
	z = 0
	traj = np.array([x,y,z])
	return traj

def IK (t): #Inverse Kinematics
	x = 0
	y = 0
	z = 0
	a1 = 22
	a2 = 125
	a3 = 150
	d1 = 71
	d2 = 0
	d3 = 0
	trajec = trajectory(t)
	x = trajec[0]
	y = trajec[1]
	z = trajec[2]

	p = z-(d1+65)
	s = math.sqrt(x*x+y*y)-30-a1
	theta1 = math.atan2(y,x) 
	theta3 = -math.acos((s*s+p*p-a2*a2-a3*a3)/(2*a2*a3))
	theta2 = math.atan2(p,s)+math.atan2(a3*math.sin(-theta3),a2+a3*math.cos(-theta3))
	theta1 = theta1*180/math.pi
	theta2 = theta2*180/math.pi
	theta3 = theta3*180/math.pi
	theta = np.array([theta1,theta2,theta3])
	return theta

def GC(pos1,pos2,pos3):	#Gravity compansation
	g1 = 0
	g2 = 937760463*math.cos(pos2+pos3)/10000 + 412323129*math.cos(pos2)
	g3 = 937760463*math.cos(pos2+pos3)/10000
	GTerms = np.array([g1,g2,g3])
	return GTerms

#-----------------
#		Init
#-----------------
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

# Motor1 Parameters
Kp1 = Kp
Ki1 = 0
Kd1 = Kd
error1 = 0
errorPrev1 = 0
errorInt1 = 0.0
errorDer1 = 0.0

# Motor2 Parameters
Kp2 = Kp
Ki2 = 0
Kd2 = Kd
error2 = 0
errorPrev2 = 0
errorInt2 = 0.0
errorDer2 = 0.0

# Motor3 Parameters
Kp3 = Kp
Ki3 = 0
Kd3 = Kd
error3 = 0
errorPrev3 = 0
errorInt3 = 0.0
errorDer3 = 0.0

DT = 0.01 # Initial step size
controlInput1 = 0
controlInput2 = 0
controlInput3 = 0

# Initial encoder readings
initPos1 = -DMCC.getQEI(0,1)*360/564
initPos2 = -DMCC.getQEI(0,2)*360/564
initPos3 = -DMCC.getQEI(1,1)*360/564

t=0 #init time
terminate=0
while True:
	DesTheta = IK(t) #Desired Theta

	while True:
		time_initial = time.time() # Record the starting time
		#--------------------		
		#	Control of Q1
		#--------------------
		pos1 = -DMCC.getQEI(0,1)*360/564 - initPos1	#Read encoder
		error1 = DesTheta[0] - pos1					#Generate the error
		errorDer1 = (error1 - errorPrev1)/DT 		#Derivative of error
		errorPrev1 = error1
		controlInput1 = Kp*error1 - Kd*errorDer1

		# Upper limit of setMotor parameter is 5000, which corresponds to 100% duty cycle for PWM signal of H-bridge, hence the saturation.
		if controlInput1 > 5000:
			controlInput1 = 5000
		elif controlInput1 < -5000:
			controlInput1 = -5000
		
		#DMCC.setMotor(cape number, motor number, power)
		DMCC.setMotor(0,1,int(controlInput1)) # Command signal to the motor

		print "Encoder-1: ", pos1, "Error-1: ", error1, "Control-1: ", (int(controlInput1)), "DT: ", DT
		#--------------------		
		#	Control of Q2
		#--------------------
		pos2 = -DMCC.getQEI(0,2)*360/564 - initPos2	#Read encoder
		error2 = DesTheta[1] - pos2					#Generate the error
		errorDer2 = (error2 - errorPrev2)/DT 		#Derivative of error
		errorPrev2 = error2
		controlInput2 = Kp*error2 - Kd*errorDer2

		# Upper limit of setMotor parameter is 5000, which corresponds to 100% duty cycle for PWM signal of H-bridge, hence the saturation.
		if controlInput2 > 5000:
			controlInput2 = 5000
		elif controlInput2 < -5000:
			controlInput2 = -5000

		#DMCC.setMotor(cape number, motor number, power)
		DMCC.setMotor(0,2,int(controlInput2)) # Command signal to the motor
		print "Encoder-2: ", pos2, "Error-2: ", error2, "Control-2: ", (int(controlInput2)), "DT: ", DT

		#--------------------		
		#	Control of Q3
		#--------------------
		pos3 = -DMCC.getQEI(1,1)*360/564 - initPos3	#Read encoder
		error3 = DesTheta[2] - pos3				#Generate the error
		errorDer3 = (error3 - errorPrev3)/DT 		#Derivative of error
		errorPrev3 = error3
		controlInput3 = Kp*error3 - Kd*errorDer3

		# Upper limit of setMotor parameter is 5000, which corresponds to 100% duty cycle for PWM signal of H-bridge, hence the saturation.
		if controlInput3 > 5000:
			controlInput3 = 5000
		elif controlInput3 < -5000:
			controlInput3 = -5000
		
		#DMCC.setMotor(cape number, motor number, power)
		DMCC.setMotor(1,1,int(controlInput3)) # Command signal to the motor
		print "Encoder-3: ", pos3, "Error-3: ", error3, "Control-3: ", (int(controlInput3)), "DT: ", DT

		#Pick the highest error
		if error1>=error2 and error1>=error3:
			terminate = error1
		if error2>=error1 and error2>=error3:
			terminate = error2
		if error3>=error2 and error3>=error2:
			terminate = error3

		if terminate < 5:
			break

		DT = time.time() - time_initial # Elapsed time in one cycle.
		t = DT + t 	#total time