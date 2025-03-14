import sys
import numpy as np
import DLM
import VLM
import copy
import json

# Step 1: Open the aerogrid JSON file
with open('aerogrid.json', 'r') as file:
    # Step 2: Load the JSON data into a Python dictionary
    aerogrid = json.load(file)
    aerogrid['N']=np.array(aerogrid['N']);
    aerogrid['l']=np.array(aerogrid['l']);
    aerogrid['A']=np.array(aerogrid['A']).T;
    aerogrid['offset_j']=np.array(aerogrid['offset_j']).T;
    aerogrid['offset_l']=np.array(aerogrid['offset_l']).T;
    aerogrid['offset_P1']=np.array(aerogrid['offset_P1']).T;
    aerogrid['offset_P3']=np.array(aerogrid['offset_P3']).T;

# Step 3: Open the user input json file
with open('user_input.json', 'r') as file:
    # Step 4: Load the JSON data into a Python dictionary
    airdata=json.load(file)
    Ma=airdata['Mach']
    k=airdata['RedFreq']
    
Qjj,Ajj,Ajj_VLM = DLM.calc_Qjj(aerogrid, Ma, k,method='quartic')
    
#Ajj is the AIC of the overall DLM solution.
#Ajj_VLM is the AIC of just the VLM part.

####wdivU_vec=np.zeros((588,1))
#####Inner flap is oscillating
####
####a=300-72-12
####for c in range(12):   #loop over spanwise strips
####    for b in range(6): # loop over chordwise strips
####        ang_amp=0.65 # in degrees 
####        distance=(0.18/6*0.75)+b*0.18/6          #because collocation points at 0.75 of chord of each box
####        wdivU_vec[a]=distance*k*np.radians(ang_amp)
####        a=a+1
####
####wdivU_vec[144:217]=wdivU_vec[216:289]
####
####cp=np.dot(np.linalg.inv(Ajj),wdivU_vec)
####L=np.dot(np.diag(aerogrid['A']),cp)
####
#####plot the correct parts of cp
####import matplotlib.pyplot as plt
####from mpl_toolkits import mplot3d
####
####
####
##### Create a 3D plot
####fig = plt.figure()
####ax = plt.axes(projection='3d')
####
##### Plot the surface
####X=aerogrid['offset_j'][:,0]
####Y=aerogrid['offset_j'][:,1]
####
####Z1=wdivU_vec.reshape(588,)
####Z2=np.real(cp.reshape(588,))
####Z3=np.real(L.reshape(588,))
####
####ax.scatter(X, Y, Z2 , 'gray')
####
##### Show the plot
####plt.show()

wdivU_vec=np.ones((1056,1))
diagS=np.diag(aerogrid['A'])
print(np.sum(np.dot(diagS,np.dot(Qjj,wdivU_vec))/np.sum(aerogrid['A'])))



        
        
