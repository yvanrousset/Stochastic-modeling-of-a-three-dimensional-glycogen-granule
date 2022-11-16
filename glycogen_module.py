#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 11:45:36 2022

@author: yvan
"""

import math
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
import os, os.path
import random
import json

class glycogen_structure:
    ''' This class generate glycogen granule. Its structural informations is stored in the information dictionary
        where keys correspond to chains name, and values are dictionary containings their identity, status, genration
        conectivity information (identity of mother and daughters) and location of the glucoses units on the chains
    '''   
    
    def __init__(self, dictio = {}, initial_number_of_chains:int = 2, l:float = 1):
        '''By default, there are two chains that point in oposite directions from the surface of a sphere
        with radius r_gn, which coresponds to the GN protein.
        l = 1, by default, is the unitless distance between two glucoses units organised in a single helix
        '''
 
        self.parameters = dictio

        r_gn = self.parameters['r_GN']
        info = {}
        for i in range(initial_number_of_chains):
            if i == 0:
                info["chain"+str(i)] = {"identity":i, "status": 2, "generation":0 , "identity_of_mother":None, "identity_of_daughter":[], "position_of_daughter":[], 
                                 "glucose_location":[[0,0,r_gn], [0,0,r_gn+1], [0,0,r_gn+2], [0,0,r_gn+3], [0,0,r_gn+4]]}
            else:
                info["chain"+str(i)] = {"identity":i, "status": 2, "generation":0 , "identity_of_mother":None,"identity_of_daughter":[], "position_of_daughter":[], 
                                 "glucose_location":[[0,0,-r_gn], [0,0,-r_gn-1] ,[0,0,-r_gn-2] ,[0,0,-r_gn-3], [0,0,-r_gn-4]]}            
        self.information = info
        self.distance_between_units = l
    

    def Get_last_index(self)->int:
        x = 0
        for chain in self.information.values():
            if chain["identity"]>0:
                x = chain["identity"] 
        return x

    def Create_a_Chain(self,index, status = None, generation = None):
        new_chain = {"identity":index, "status": status, "generation":generation , 
                    "identity_of_mother":None, "identity_of_daughter":[], "position_of_daughter":[], 
                    "glucose_location":[]}
        #print(index)
        self.information['chain'+str(index)] = new_chain
    
    def write_structure(self, name = "glycogen_structure.json"):
        with open(name, "w") as f:
            json.dump(self.information, f, indent=2)

    def number_of_glucose_fixed(self):
        compteur = 0
        for chain in self.information.values():
            compteur= compteur + len(chain['glucose_location'])
        return compteur

    def number_of_chain(self):
        return len(self.information)

    ###############################################
    # All methods that "find" something in the data
    def Find_chain_for_gs(self) -> list:
        '''This function look into the information of the structure and return the chains that are possible substrate for gs
        It is done by returning the identity of the chains that dont have status '0' which correspond to dp that are to short to bind GS
        '''
        list_chains_substrate_for_gs = [];        
        for chain in self.information.values():
            if chain["status"] != 0:
                list_chains_substrate_for_gs.append(chain["identity"])
        return list_chains_substrate_for_gs
    

    def Find_chain_for_gbe(self) -> list:
        list_chains_substrate_for_gbe = [];
        for chain in self.information.values():
            if chain["status"] == 3:
                list_chains_substrate_for_gbe.append(chain["identity"])
        return list_chains_substrate_for_gbe 

    def Find_chain_for_gp(self) -> list:
        list_chains_substrate_for_gp = [];
        for chain in self.information.values():
            if chain["status"] in [2, 3]:
               list_chains_substrate_for_gp.append(chain["identity"])
        return list_chains_substrate_for_gp
    
    def Find_chain_for_gde(self) -> list:
        list_chains_substrate_for_gde = []     
        for chain in self.information.values():
            if chain["status"] == 1:
                list_chains_substrate_for_gde.append(chain["identity"])           
        return list_chains_substrate_for_gde
    #########################################################           
    
    def Get_chain_from_identity(self,identity):
        for chain in self.information.values(): 
            if chain["identity"] == identity:
                return chain
        raise AttributeError('no chain matchs identity')
        
                
                
    def Scan_excluded_volume(self, position_to_test, identity=None):
        #print((position_to_test[0]-self.parameters["r_GN"])**2 + (position_to_test[1])**2 + (position_to_test[2])**2,self.parameters["r_GN"]**2)
        #print("here is the position to test:",position_to_test)
        if (position_to_test[0])**2 + (position_to_test[1])**2 + (position_to_test[2])**2 < self.parameters["r_GN"]**2:
            return False        
        
        for chain in self.information.values():  
            glucose_position = chain["glucose_location"]
            for pos in glucose_position:
                if (position_to_test[0]-pos[0])**2 + (position_to_test[1]-pos[1])**2 + (position_to_test[2]-pos[2])**2 < self.parameters["b"]**2:
                    if identity != chain["identity"]:
                        return False
        return True
    
    def Torsion_function(self):
        pass

    def Gillespie_step_ma(self) -> tuple:
        ''' This functions takes concentrations of the enzymes and the structure info of a glycogen granules and
        return what is the next reaction to occurs and which time has been spent. (Following a gillespie algorithm)
        Assuming mass action kinetics
        '''

        #propensity assuming mass action kinetics
        h_gs  = self.parameters['GS']*len(self.Find_chain_for_gs())
        h_gp  = self.parameters["GP"]*len(self.Find_chain_for_gp())
        h_gbe = self.parameters["GBE"]*len(self.Find_chain_for_gbe())
        h_gde = self.parameters["GDE"]*len(self.Find_chain_for_gde())
        
        a = h_gs + h_gp + h_gbe + h_gde
        
        if a == 0:
            return "no reaction can be proceed, all propensities are zero",0
        r2=random.uniform(0,a)    
        r1=random.uniform(0,1)	
        
        d_t = (1/a)*math.log(1/r1)

        if r2 < h_gs :
            return ("Act_gs",d_t)
        if r2 >= h_gs and r2 < h_gs + h_gp :
            return ("Act_gp",d_t)
        if r2 >=  h_gs + h_gp  and r2 < h_gs + h_gp + h_gbe :
            return ("Act_gbe",d_t)
        if r2 >=  h_gs + h_gp + h_gbe  and r2 < h_gs + h_gp + h_gbe + +h_gde:
            return ("Act_gde",d_t)          


    ##############################################
    # Functions that try to perform a change in the structure after a reactions was picked
    def Act_gs(self, selected_chain_index = None)-> bool:
        if selected_chain_index == None:
            selected_chain_index = random.choice(self.Find_chain_for_gs())
                
        chain = self.Get_chain_from_identity(selected_chain_index)
        
        glucose_position = chain["glucose_location"]
        
        theta, phi = angle3Dchain(glucose_position)           
        l = self.distance_between_units                                         
        x_nre, y_nre, z_nre = glucose_position[-1]

        x=x_nre+l*math.cos(phi)*math.cos(theta);                   # new x monomer position    
        y=y_nre+l*math.cos(phi)*math.sin(theta);                   # new y monomer position
        z=z_nre+l*math.sin(phi)
        
        if self.Scan_excluded_volume([x,y,z],identity = chain["identity"]) is True:
            self.information["chain"+str(selected_chain_index)]['glucose_location'].append([x,y,z])
            self.update_chain_status(chain["identity"])
            return True
        else:
            return False    
                
    def Act_gp(self)-> bool:
        try:
            selected_chain_index = random.choice(self.Find_chain_for_gp())    
            #print(selected_chain_index)        
            chain = self.Get_chain_from_identity(selected_chain_index)        
            del self.information["chain"+str(selected_chain_index)]['glucose_location'][-1] 
            self.update_chain_status(chain['identity'])
            return True
        
        except IndexError:
            return False 
  
    
    def Act_gbe_flexible_model(self)-> bool:  
        if self.Find_chain_for_gbe() == []:
            return False
        
        else:
            selected_chain_index = random.choice(self.Find_chain_for_gbe())    
            chain = self.Get_chain_from_identity(selected_chain_index) 
            eps = 5.4

            #selection of the position to cleave:
            N = len(chain["glucose_location"])

            if len(chain["position_of_daughter"])==0:
                ind = 0
            else:
                ind = chain["position_of_daughter"][-1]+1
                
            sub_list_to_cleave = [k for k in range(
                                ind+self.parameters['size_spec_gbe_spacing'] 
                                +self.parameters['size_spec_gbe_leftover'],
                                N-self.parameters['size_spec_gbe_transferred'])]  

            #if sub_list_to_cleave !=0:
            indice_where_to_cleave = random.choice(sub_list_to_cleave)
            #print(sub_list_to_cleave)
            N_new = N-indice_where_to_cleave-1
            #print('new',N_new)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            
            #print(sub_list_to_branch)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,ind+1+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'])]  
            #if sub_list_to_branch !=0:
            indice_where_to_branch = random.choice(sub_list_to_branch)
            
            #print(indice_where_to_cleave, indice_where_to_branch)

            chain_length = N-(indice_where_to_cleave+1)

            alpha0, beta0 = angle3Dchain( chain["glucose_location"] ) 

            theta = random.uniform(-pi,pi)
            phi = random.uniform(-pi/2,pi/2)


            x0,y0,z0 = chain["glucose_location"][indice_where_to_branch]

            l = self.distance_between_units

            xstart=x0+(l+eps)*math.cos(phi+beta0)*math.cos(theta+alpha0);                                   #position on x axis of the first monomer on the new chain
            ystart=y0+(l+eps)*math.cos(phi+beta0)*math.sin(theta+alpha0);
            zstart=z0+(l+eps)*math.sin(phi+beta0)

            xend=xstart+(chain_length-1)*l*math.cos(phi+beta0)*math.cos(theta+alpha0);                              #position on x axis of the last monomer on the new chain
            yend=ystart+(chain_length-1)*l*math.cos(phi+beta0)*math.sin(theta+alpha0);                              #position on y axis of the last monomer on the new chain
            zend=zstart+(chain_length-1)*l*math.sin(phi+beta0)

            X=np.linspace(xstart,xend,chain_length) 
            Y=np.linspace(ystart,yend,chain_length)
            Z=np.linspace(zstart,zend,chain_length)

            cond = 0
            for k in range(chain_length):
                if self.Scan_excluded_volume([X[k],Y[k],Z[k]]) == False:
                    cond = 1
            if cond == 0:
                #creation of the new chain:
                index = self.Get_last_index()
                self.Create_a_Chain(index+1)

                self.information['chain'+str(index+1)]['generation'] = chain["generation"]+1
                for k in range(chain_length):
                    self.information['chain'+str(index+1)]['glucose_location'].append([X[k],Y[k],Z[k]]) 
                self.information['chain'+str(index+1)]['identity_of_mother'] = selected_chain_index
                new_chain = self.information['chain'+str(index+1)]

                self.information['chain'+str(selected_chain_index)]['identity_of_daughter'].append(index+1)
                self.information['chain'+str(selected_chain_index)]['position_of_daughter'].append(indice_where_to_branch)
                del self.information['chain'+str(selected_chain_index)]['glucose_location'][indice_where_to_cleave+1:] 
                
                self.update_chain_status(new_chain['identity'])
                self.update_chain_status(chain['identity'])
                #print(self.information)
                return True

            else :
                return False             



    def Act_gbe_flexible_model_new(self)-> bool:  
        '''This function try to perform a branching on the glycogen structure. It returns True if the branching has been made and False if not (due to steric hindrance conflict).
        The branching geometry differs from Act_gbe(). 
        '''

        if self.Find_chain_for_gbe() == []:
            return False
        
        else:
            selected_chain_index = random.choice(self.Find_chain_for_gbe())    
            chain = self.Get_chain_from_identity(selected_chain_index) 
            eps = 5.4

            #selection of the position to cleave:
            N = len(chain["glucose_location"])

            if len(chain["position_of_daughter"])==0:
                ind = 0
            else:
                ind = chain["position_of_daughter"][-1]+1
                
            sub_list_to_cleave = [k for k in range(
                                ind+self.parameters['size_spec_gbe_spacing'] 
                                +self.parameters['size_spec_gbe_leftover'],
                                N-self.parameters['size_spec_gbe_transferred'])] 
            #if sub_list_to_cleave !=0:
            indice_where_to_cleave = random.choice(sub_list_to_cleave)
            #print(sub_list_to_cleave)
            N_new = N-indice_where_to_cleave-1
            #print('new',N_new)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            
            #print(sub_list_to_branch)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,ind+1+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'])]  
            #if sub_list_to_branch !=0:
            indice_where_to_branch = random.choice(sub_list_to_branch)
            
            #print(indice_where_to_cleave, indice_where_to_branch)

            chain_length = N-(indice_where_to_cleave+1)
            #alpha0, beta0 = angle3Dchain( chain["glucose_location"] ) 

            dir_x, dir_y, dir_z = chain_normed_vector( chain["glucose_location"] ) 
            [xstart,ystart,zstart] = generate_perpendicular_vector([dir_x, dir_y, dir_z ])
            #print('direction=',(dir_x, dir_y, dir_z),'scalar product= ',np.dot([dir_x, dir_y, dir_z],[xstart,ystart,zstart]))

            x0,y0,z0 = chain["glucose_location"][indice_where_to_branch]
            l = self.distance_between_units

            xstart=x0+(l+eps)*xstart                                  #position on x axis of the first monomer on the new chain
            ystart=y0+(l+eps)*ystart
            zstart=z0+(l+eps)*zstart

            theta = random.uniform(-math.pi,math.pi)
            delta = random.uniform(-math.pi/2,math.pi/2)

            xend=xstart+(chain_length-1)*l*math.cos(delta)*math.cos(theta);                              #position on x axis of the last monomer on the new chain
            yend=ystart+(chain_length-1)*l*math.cos(delta)*math.sin(theta);                              #position on y axis of the last monomer on the new chain
            zend=zstart+(chain_length-1)*l*math.sin(delta)

            X=np.linspace(xstart,xend,chain_length) 
            Y=np.linspace(ystart,yend,chain_length)
            Z=np.linspace(zstart,zend,chain_length)

            cond = 0
            for k in range(chain_length):
                if self.Scan_excluded_volume([X[k],Y[k],Z[k]]) == False:
                    cond = 1
            if cond == 0:
                #creation of the new chain:
                index = self.Get_last_index()
                self.Create_a_Chain(index+1)

                self.information['chain'+str(index+1)]['generation'] = chain["generation"]+1
                for k in range(chain_length):
                    self.information['chain'+str(index+1)]['glucose_location'].append([X[k],Y[k],Z[k]]) 
                self.information['chain'+str(index+1)]['identity_of_mother'] = selected_chain_index
                new_chain = self.information['chain'+str(index+1)]
                #print(new_chain)

                #print('identities',(new_chain['identity'],chain['identity']))

            #update previous chain

                self.information['chain'+str(selected_chain_index)]['identity_of_daughter'].append(index+1)
                self.information['chain'+str(selected_chain_index)]['position_of_daughter'].append(indice_where_to_branch)
                del self.information['chain'+str(selected_chain_index)]['glucose_location'][indice_where_to_cleave+1:] 
                
                self.update_chain_status(new_chain['identity'])
                self.update_chain_status(chain['identity'])
                #print(self.information)
                return True

            else :
                return False 

    def Act_gbe_strict_location(self)-> bool:  
        if self.Find_chain_for_gbe() == []:
            return False
        
        else:
            selected_chain_index = random.choice(self.Find_chain_for_gbe())    
            chain = self.Get_chain_from_identity(selected_chain_index) 
            eps = 5.4

            #selection of the position to cleave:
            N = len(chain["glucose_location"])

            if len(chain["position_of_daughter"])==0:
                ind = 0
            else:
                ind = chain["position_of_daughter"][-1]+1
                
            sub_list_to_cleave = [k for k in range(
                                ind+self.parameters['size_spec_gbe_spacing'] 
                                +self.parameters['size_spec_gbe_leftover'],
                                N-self.parameters['size_spec_gbe_transferred'])]  

            #if sub_list_to_cleave !=0:
            indice_where_to_cleave = sub_list_to_cleave[-1]
            #print(sub_list_to_cleave)
            N_new = N-indice_where_to_cleave-1
            #print('new',N_new)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            
            #print(sub_list_to_branch)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,ind+1+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'])]  
            #if sub_list_to_branch !=0:
            indice_where_to_branch = sub_list_to_branch[-1]
            
            #print(indice_where_to_cleave, indice_where_to_branch)

            chain_length = N-(indice_where_to_cleave+1)

            alpha0, beta0 = angle3Dchain( chain["glucose_location"] ) 

            theta = random.uniform(-pi,pi)
            phi = random.uniform(-pi/2,pi/2)


            x0,y0,z0 = chain["glucose_location"][indice_where_to_branch]

            l = self.distance_between_units

            xstart=x0+(l+eps)*math.cos(phi+beta0)*math.cos(theta+alpha0);                                   #position on x axis of the first monomer on the new chain
            ystart=y0+(l+eps)*math.cos(phi+beta0)*math.sin(theta+alpha0);
            zstart=z0+(l+eps)*math.sin(phi+beta0)

            xend=xstart+(chain_length-1)*l*math.cos(phi+beta0)*math.cos(theta+alpha0);                              #position on x axis of the last monomer on the new chain
            yend=ystart+(chain_length-1)*l*math.cos(phi+beta0)*math.sin(theta+alpha0);                              #position on y axis of the last monomer on the new chain
            zend=zstart+(chain_length-1)*l*math.sin(phi+beta0)

            X=np.linspace(xstart,xend,chain_length) 
            Y=np.linspace(ystart,yend,chain_length)
            Z=np.linspace(zstart,zend,chain_length)

            cond = 0
            for k in range(chain_length):
                if self.Scan_excluded_volume([X[k],Y[k],Z[k]]) == False:
                    cond = 1
            if cond == 0:
                #creation of the new chain:
                index = self.Get_last_index()
                self.Create_a_Chain(index+1)

                self.information['chain'+str(index+1)]['generation'] = chain["generation"]+1
                for k in range(chain_length):
                    self.information['chain'+str(index+1)]['glucose_location'].append([X[k],Y[k],Z[k]]) 
                self.information['chain'+str(index+1)]['identity_of_mother'] = selected_chain_index
                new_chain = self.information['chain'+str(index+1)]
                #print(new_chain)

                #print('identities',(new_chain['identity'],chain['identity']))

            #update previous chain

                self.information['chain'+str(selected_chain_index)]['identity_of_daughter'].append(index+1)
                self.information['chain'+str(selected_chain_index)]['position_of_daughter'].append(indice_where_to_branch)
                del self.information['chain'+str(selected_chain_index)]['glucose_location'][indice_where_to_cleave+1:] 
                
                self.update_chain_status(new_chain['identity'])
                self.update_chain_status(chain['identity'])
                #print(self.information)
                return True

            else :
                return False       

    def Act_gbe_strict_location_variation(self)-> bool:  
        if self.Find_chain_for_gbe() == []:
            return False
        
        else:
            selected_chain_index = random.choice(self.Find_chain_for_gbe())    
            chain = self.Get_chain_from_identity(selected_chain_index) 
            eps = 5.4

            #selection of the position to cleave:
            N = len(chain["glucose_location"])

            if len(chain["position_of_daughter"])==0:
                ind = 0
            else:
                ind = chain["position_of_daughter"][-1]+1
                
            sub_list_to_cleave = [k for k in range(
                                ind+self.parameters['size_spec_gbe_spacing'] 
                                +self.parameters['size_spec_gbe_leftover'],
                                N-self.parameters['size_spec_gbe_transferred'])]  

            #if sub_list_to_cleave !=0:
            indice_where_to_cleave = sub_list_to_cleave[0]
            #print(sub_list_to_cleave)
            N_new = N-indice_where_to_cleave-1
            #print('new',N_new)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,N-N_new-self.parameters['size_spec_gbe_leftover'])]  
            
            #print(sub_list_to_branch)
            #sub_list_to_branch = [k for k in range(ind+self.parameters['size_spec_gbe_spacing'] ,ind+1+self.parameters['size_spec_gbe_spacing']+self.parameters['size_spec_gbe_leftover'])]  
            #if sub_list_to_branch !=0:
            indice_where_to_branch = sub_list_to_branch[-1]
            
            #print(indice_where_to_cleave, indice_where_to_branch)

            chain_length = N-(indice_where_to_cleave+1)

            alpha0, beta0 = angle3Dchain( chain["glucose_location"] ) 

            theta = random.uniform(-pi,pi)
            phi = random.uniform(-pi/2,pi/2)


            x0,y0,z0 = chain["glucose_location"][indice_where_to_branch]

            l = self.distance_between_units

            xstart=x0+(l+eps)*math.cos(phi+beta0)*math.cos(theta+alpha0);                                   #position on x axis of the first monomer on the new chain
            ystart=y0+(l+eps)*math.cos(phi+beta0)*math.sin(theta+alpha0);
            zstart=z0+(l+eps)*math.sin(phi+beta0)

            xend=xstart+(chain_length-1)*l*math.cos(phi+beta0)*math.cos(theta+alpha0);                              #position on x axis of the last monomer on the new chain
            yend=ystart+(chain_length-1)*l*math.cos(phi+beta0)*math.sin(theta+alpha0);                              #position on y axis of the last monomer on the new chain
            zend=zstart+(chain_length-1)*l*math.sin(phi+beta0)

            X=np.linspace(xstart,xend,chain_length) 
            Y=np.linspace(ystart,yend,chain_length)
            Z=np.linspace(zstart,zend,chain_length)

            cond = 0
            for k in range(chain_length):
                if self.Scan_excluded_volume([X[k],Y[k],Z[k]]) == False:
                    cond = 1
            if cond == 0:
                #creation of the new chain:
                index = self.Get_last_index()
                self.Create_a_Chain(index+1)

                self.information['chain'+str(index+1)]['generation'] = chain["generation"]+1
                for k in range(chain_length):
                    self.information['chain'+str(index+1)]['glucose_location'].append([X[k],Y[k],Z[k]]) 
                self.information['chain'+str(index+1)]['identity_of_mother'] = selected_chain_index
                new_chain = self.information['chain'+str(index+1)]
                #print(new_chain)

                #print('identities',(new_chain['identity'],chain['identity']))

            #update previous chain

                self.information['chain'+str(selected_chain_index)]['identity_of_daughter'].append(index+1)
                self.information['chain'+str(selected_chain_index)]['position_of_daughter'].append(indice_where_to_branch)
                del self.information['chain'+str(selected_chain_index)]['glucose_location'][indice_where_to_cleave+1:] 
                
                self.update_chain_status(new_chain['identity'])
                self.update_chain_status(chain['identity'])
                #print(self.information)
                return True

            else :
                return False  


    def Act_gde(self)-> bool:
        if self.Find_chain_for_gde() == []:
            return False
        selected_chain_index = random.choice(self.Find_chain_for_gde()) 
        chain = self.Get_chain_from_identity(selected_chain_index)
        N = len(chain['glucose_location'])

        for i in range(N-1):
            self.Act_gs(chain['identity_of_mother'])
        
        del self.information['chain'+str(selected_chain_index)]
        del self.information['chain'+str(chain['identity_of_mother'])]['identity_of_daughter'][-1]
        del self.information['chain'+str(chain['identity_of_mother'])]['position_of_daughter'][-1]
        return True
        
    ##############################################

    def update_chain_status(self, identity):
        #print(chain["glucose_location"])
        
        chain = self.Get_chain_from_identity(identity)
        N = len(chain["glucose_location"])
        if len(chain["position_of_daughter"])==0:
            ind = 0
        else:
            ind = chain["position_of_daughter"][-1]+1
            
        if N-ind < self.parameters['size_spec_gys'] :
            self.information['chain'+str(identity)]["status"] = 1
        elif N-ind == self.parameters['size_spec_gys'] and len(self.information['chain'+str(identity)]["identity_of_daughter"])==0 :
            self.information['chain'+str(identity)]["status"] = 1
        elif N-ind <= self.parameters['size_spec_gbe_spacing'] + self.parameters['size_spec_gbe_leftover'] + self.parameters['size_spec_gbe_transferred'] :
            self.information['chain'+str(identity)]["status"] = 2
            
        else:
            self.information['chain'+str(identity)]["status"] = 3
            

    def plot_structure(self):
        XLIST=[];YLIST=[];ZLIST=[]
        Xnre=[];Ynre=[];Znre=[]
        Xa = []
        Ya = []
        Za = []
        alpha_segment = []
        for chain in self.information.values():
            for pos in chain['glucose_location']: 
                XLIST.append(pos[0])
                YLIST.append(pos[1])
                ZLIST.append(pos[2])

            for branching_index, daugther_index in zip(chain['position_of_daughter'],chain['identity_of_daughter']):
                #print('identity:',chain['identity'],branching_index, daugther_index)
                
                x1 = chain['glucose_location'][branching_index]
                daughter_chain = self.Get_chain_from_identity(daugther_index)    
                #print('position_of_daughter',daughter_chain['position_of_daughter'])              
                x2 = daughter_chain['glucose_location'][0]

                alpha_segment.append([x1,x2])

                # if chain['identity'] == 0:
                #     Xa.append(chain['glucose_location'][2][0])
                #     Ya.append(chain['glucose_location'][2][1])
                #     Za.append(chain['glucose_location'][2][2])
                # else:
                #     Xa.append(chain['glucose_location'][0][0])
                #     Ya.append(chain['glucose_location'][0][1])
                #     Za.append(chain['glucose_location'][0][2])


        XLIST,YLIST,ZLIST= 0.24*np.asarray(XLIST),  0.24*np.asarray(YLIST), 0.24* np.asarray(ZLIST)
        ax=plt.figure().gca(projection='3d')
        ax.plot(XLIST,YLIST,ZLIST,'o',markersize=7, markeredgewidth=0.2,markeredgecolor='black',label='glycogen 3d')
        #ax.view_init(elev=45, azim=45)
        ax.view_init(elev=15, azim=60)
        plt.xlabel('x [nm]',fontsize='15')
        #ax.plot(0.24*np.asarray(Xa),0.24*np.asarray(Ya),0.24*np.asarray(Za),'-',markersize=1, color ='red')
        #print(alpha_segment)
        plt.ylabel('y [nm]',fontsize='15')
        ax.set(xlim=(-30, 30), ylim=(-30, 30),zlim=(-30,30))
        #plt.savefig('test.elongated.png',quality=95)

        #for i in range(alpha_segment):
        #    ax.plot(0.24*np.asarray([alpha_segment[i][0][0],alpha_segment[i][1][0]]),0.24*np.asarray([alpha_segment[i][0][1],alpha_segment[i][1][1]]),0.24*np.asarray([alpha_segment[i][0][2],alpha_segment[i][1][2]]),'-',markersize=1, color ='red')
        
        plt.show()
        ax.legend()


    def get_radius(self, unit = 'adim'):
        listex=[];listey=[];listez=[]

        if unit == 'adim':
            alpha = 1.0
        else:
            alpha = 0.24
        for chain in self.information.values():  
            glucose_position = chain["glucose_location"]
            for pos in glucose_position:
                listex.append(pos[0])
                listey.append(pos[1])
                listez.append(pos[2])
        
        N  = len(listex)
        xmean=1./N*sum(listex)
        ymean=1./N*sum(listey)
        zmean=1./N*sum(listez)
        
        r2=1.0/N*sum((listex-xmean*np.ones(N))**2+(listey-ymean*np.ones(N))**2+(listez-zmean*np.ones(N))**2)
        return alpha*r2**0.5, [xmean,ymean,zmean]       


    def AtoBratio(self):
        A = 0
        B = 0
        for chain in self.information.values():  
            if len(chain["identity_of_daughter"]) == 0:
                A += 1
            else:
                B += 1
        return A/B

    def avg_cl(self):
        cl_list = []
        for chain in self.information.values(): 
            cl_list.append(len(chain["glucose_location"]))
        return np.mean(cl_list)       


    def get_last_gen_index(self):
        max_gen = 0
        for chain in self.information.values():
            gen = chain['generation']
            if gen > max_gen:
                max_gen = gen
        return max_gen

    def bd(self):
        number_16 = self.number_of_chain()-1
        number_14 =  self.number_of_chain()*(self.avg_cl()-1)
        return number_16/number_14

    def occupancy(self):
        v = 1.33*0.24
        r_g,_ = self.get_radius(unit = 'nm')
        return self.number_of_glucose_fixed()*v/(4/3*math.pi*r_g**3)


    def get_number_of_units_within_radius(self,radius = 10):
        _, [x,y,z] = self.get_radius()
        counter = 0
        for chain in self.information.values():
            for pos in chain['glucose_location']:
                if ((x-pos[0])**2+(y-pos[1])**2+(z-pos[2])**2)**0.5< radius:
                    counter+=1
        return counter 


    def get_density_profile(self, r_sample = np.linspace(1,100,50), delta_r = 5):
        r_tot, center = self.get_radius()
        occupancy_list = []

        for r in r_sample:
            N_ext = self.get_number_of_units_within_radius(r+delta_r)
            N_in  = self.get_number_of_units_within_radius(r)

            Volume_layer_in_nm3 = (4.0/3) * math.pi * (r*0.24 + delta_r*0.24)**3 - (4.0/3) * math.pi * (r*0.24)**3
            occupancy = (1.33*0.24)*(N_ext-N_in)/Volume_layer_in_nm3

            occupancy_list.append(occupancy)

        return occupancy_list


    def time_course(self, time, external_parameters):
        pass



# Some other usefull functions
def chain_normed_vector(liste):
    return np.asarray(liste[-1])-np.asarray(liste[-2])

def angle3Dchain(liste):    
    X0=np.asarray(liste[-2])
    X1=np.asarray(liste[-1])
    [x,y,z]=X1-X0;

    length = (x*x + y*y + z*z)**0.5
    phi=math.asin(z/length)
    
    if y>=0 and (x*x+y*y!=0):
        theta=math.acos(x/(x*x+y*y)**0.5)
    elif y>=0 and (x*x+y*y==0):
        theta=0.0;
    else:
        theta=2*math.pi-math.acos(x/(x*x+y*y)**0.5)        
    return theta,phi       



def generate_perpendicular_vector(u):
    '''This function take a 3d vector u, and generate a random normed perpendicular vector, v to u
    '''
    theta = random.uniform(-math.pi,math.pi)
    delta = random.uniform(-math.pi/2,math.pi/2)

    [xrand,yrand,zrand] = [math.cos(delta)*math.cos(theta), math.cos(delta)*math.sin(theta), math.sin(delta)]
    
    w = np.array([xrand, yrand, zrand])
    u = np.array(u)

    a = np.dot(w,u) / np.dot(u,u)
    v = w-np.dot(a,u)
    normv = np.linalg.norm(v)

    return(v/normv)








""" 
            #creation of the new chain:
            index = self.Get_last_index()
            self.Create_a_Chain(index+1)

            self.information['chain'+str(index+1)]['generation'] = chain["generation"]+1
            for k in range(chain_length):
                self.information['chain'+str(index+1)]['glucose_location'].append([X[k],Y[k],Z[k]]) 
            self.information['chain'+str(index+1)]['identity_of_mother'] = selected_chain_index
            new_chain = self.information['chain'+str(index+1)]

            #update previous chain

            self.information['chain'+str(selected_chain_index)]['identity_of_daughter'].append(index+1)
            self.information['chain'+str(selected_chain_index)]['position_of_daughter'].append(indice_where_to_branch)
            del self.information['chain'+str(selected_chain_index)]['glucose_location'][indice_where_to_cleave+1:] 

            #print('identities',(new_chain['identity'],chain['identity']))
            self.update_chain_status(new_chain['identity'])
            self.update_chain_status(chain['identity'])

            return True """