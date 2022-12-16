import logging
from glycogen_module.core import GlycogenStructure
import numpy as np

from pylab import pi, plt
logger = logging.getLogger(__name__)

def plot_structure(gs: GlycogenStructure):
        XLIST=[];YLIST=[];ZLIST=[]

        for chain in gs.chains:
            for pos in chain.glucose_positions:
                #logger.debug(f"chain {chain.id}, glucose_positions: {chain.glucose_positions}.Pos: {pos}") 
                XLIST.append(pos[0])
                YLIST.append(pos[1])
                ZLIST.append(pos[2])


        XLIST,YLIST,ZLIST= GlycogenStructure.L*np.asarray(XLIST),  GlycogenStructure.L*np.asarray(YLIST), GlycogenStructure.L* np.asarray(ZLIST)
        #ax=plt.figure().gca(projection='3d')

        ax=plt.figure()
        ax = ax.add_subplot(projection='3d')


        ax.plot(XLIST,YLIST,ZLIST,'o',markersize=7, markeredgewidth=0.2,markeredgecolor='black',label='glycogen 3d')
        ax.view_init(elev=15, azim=60)
        plt.xlabel('x [nm]',fontsize='15')

        plt.ylabel('y [nm]',fontsize='15')
        ax.set(xlim=(-30, 30), ylim=(-30, 30),zlim=(-30,30))

        plt.show()
        ax.legend()