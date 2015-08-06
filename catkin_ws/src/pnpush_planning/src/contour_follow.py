#!/usr/bin/env python

# Peter KT Yu, Aug 2015
import sys
import numpy as np
from ik.roshelper import ROS_Wait_For_Msg
from geometry_msgs.msg import WrenchStamped
import tf
import tf.transformations as tfm

def main():
    setCart = rospy.ServiceProxy('/robot2_SetCartesian', robot_SetCartesian)
    
    listener = tf.TransformListener()
    
    z = 0.30
    limit = 10000
    ori = [0, 0.7071, 0.7071, 0]
    threshold = 0.01  # need to be tuned
    probe_rad = 0.00626/2
    
    # 0. move to startPos
    start_pos = [0, 0.4, z + 0.1]
    setCart(*start_pos, *ori)
    
    start_pos = [0, 0.4, z]
    setCart(*start_pos, *ori)
    curr_pos = start_pos
    
    # 1. move in -y direction till contact -> normal
    direc = np.array([0,-0.001,0])
    normal = [0,0,0]
    while True:
        curr_pos = np.array(curr_pos) + direc
        setCart(curr_pos)
        ft = ROS_Wait_For_Msg('/netft_data', sensor_msgs.msg.JointState).getmsg()
        
        (pos,ori) = listener.lookupTransform('base_link','link_ft',t)
        
        rotmat = tfm.quaternion_matrix(ori)
        w = ft.wrench
        ft_global = np.dot(rotmat, [w.x, w.y, w.z, 1.0])
        
        # in contact
        if np.sqrt(ft_global[0]**2 + ft_global[1]**2) > threshold:
            normal = ft_global[0:3] / sqrt(dot(ft_global, ft_global))
            break
        
    
    # 2. use the normal to move along the block
    index = 0
    contact_pt = []
    contact_nm = []
    while True:
        # 2.1 move 
        direc = dot(tfm.euler_matrix(0,0,1.57) , normal)[0:3]
        curr_pos = np.array(curr_pos) + direc
        setCart(curr_pos)
        # 
        (pos,ori) = listener.lookupTransform('base_link','link_ft',t)
        rotmat = tfm.quaternion_matrix(ori)
        w = ft.wrench
        
        ft_global = np.dot(rotmat, [w.x, w.y, w.z, 1.0])
        if np.sqrt(ft_global[0]**2 + ft_global[1]**2) > threshold:
            # ft_global[2] = 0  ??
            normal = ft_global[0:3] / sqrt(dot(ft_global, ft_global))
            contact_nm.append(normal)
            contact_pt.append(curr_pos + normal * probe_rad)
        
        if len(contact_pt) >= limit:
            break
    
    
if __name__=='__main__':
    main(sys.argv)