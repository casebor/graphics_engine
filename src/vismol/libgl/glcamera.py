#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import vismol.utils.matrix_operations as mop
from logging import getLogger

logger = getLogger(__name__)


class GLCamera():
    """ The GLCamera object creates a "camera" to be used in OpenGL.
        It automatically creates a viewing and projection matrices with the
        values defined in the constructor.
    """
    
    def __init__ (self, fov=30.0, var=(4.0 / 3.0), pos=np.array([0,0,10], dtype=np.float32),
                  zrp=np.array([0,0,0], dtype=np.float32)):
        """ Depending of the distance from the camera position to a defined
            reference point, it creates different clipping planes (this function
            could be improved, but will work for now).
            
            Input parameters:
                fov -- Specifies the field of view angle, in degrees, in the
                       Y direction.
                var -- Specifies the aspect ratio that determines the field of
                       view in the x direction. The aspect ratio is the ratio
                       of x (width) to y (height).
                pos -- The position in world coordinates of the camera.
                zrp -- The zero reference point defined for rotation functions.
            
            Automatically generates:
                z_near -- Specifies the distance from the viewer to the near
                          clipping plane (always positive).
                z_far -- Specifies the distance from the viewer to the far
                         clipping plane (always positive).
                fog_start -- Specifies the beggining of the fog effect (always
                             positive and lower than z_far).
                fog_end -- Specifies the end of the fog effect (always positive
                           and preferable equal to z_far).
                view_matrix -- The viewing matrix used in the render by shaders.
                projection_matrix -- The projection matrix used in the render
                                     by shaders.
        """
        self.field_of_view = np.float32(fov)
        self.viewport_aspect_ratio = np.float32(var)
        self.zero_reference_point = np.array(zrp, dtype=np.float32)
        self.max_vertical_angle = 85.0  # must be less than 90 to avoid gimbal lock
        self.horizontal_angle = 0.0
        self.vertical_angle = 0.0
        self.min_znear = 0.1
        self.min_zfar = 9.0
        dist = np.linalg.norm(pos - zrp)
        if dist <= 10.0:
            self.z_near = dist - 3.0
            self.z_far = dist + 3.0
        elif dist <= 20.0:
            self.z_near = dist - 6.0
            self.z_far = dist + 6.0
        elif dist <= 40.0:
            self.z_near = dist - 9.0
            self.z_far = dist + 9.0
        elif dist <= 80.0:
            self.z_near = dist - 12.0
            self.z_far = dist + 12.0
        else:
            self.z_near = dist - 15.0
            self.z_far = dist + 15.0
        
        self.fog_end = self.z_far
        self.fog_start = self.fog_end - self.min_zfar
        self.view_matrix = self._get_view_matrix(pos)
        self.projection_matrix = self._get_projection_matrix()
    
    def _get_view_matrix(self, position):
        """ Creates the view matrix, i.e. the matrix with the position and
            orientation of the camera used in OpenGL.
            First creates a translation matrix at the position defined by the
            "position" argument, then creates a orientation matrix initially
            with the vertical angle, then the horizontal angle is applied to
            this matrix. Finally the multiplication of translation x orientation
            is returned as the view matrix.
            
            Input parameters:
                position -- a numpy array of size 3 and ndim of 1 containing
                            the camera's position in XYZ coordinates.
            
            Returns:
                view -- a numpy array of size 16 and ndim of 2 containing the
                        information for the camera's position and orientation.
        """
        trans = mop.my_glTranslatef(np.identity(4,dtype=np.float32), -position)
        orient = mop.my_glRotatef(np.identity(4,dtype=np.float32), self.vertical_angle, np.array([1,0,0]))
        orient = mop.my_glRotatef(orient, self.horizontal_angle, np.array([0,1,0]))
        view = mop.my_glMultiplyMatricesf(trans, orient)
        return view
    
    def _get_projection_matrix(self):
        """ Creates the projection matrix using the parameters supplied in the
            constructor of the class. This function creates a perspective with
            a defined field of view, viewport aspect ratio, near and far depth
            clipping planes.
            
            Returns:
                persp -- a numpy array of size 16 and ndim of 2, containing the
                         information for the camera's field and depth view.
        """
        assert(self.field_of_view>0.0 and self.field_of_view<180.0)
        assert(self.z_near>0.0)
        assert(self.z_near<self.z_far)
        assert(self.viewport_aspect_ratio>0.0)
        persp = mop.my_glPerspectivef(self.field_of_view,self.viewport_aspect_ratio,self.z_near,self.z_far)
        return persp
    
    def get_position(self):
        """ Returns the x, y, z position of the camera in 
            absolute coordinates.
        """
        return mop.get_xyz_coords(self.view_matrix)
    
    def get_modelview_position(self, model_matrix):
        modelview = mop.my_glMultiplyMatricesf(model_matrix, self.view_matrix)
        crd_xyz = -1 * np.asmatrix(modelview[:3,:3]) * np.asmatrix(modelview[3,:3]).T
        return crd_xyz.A1
    
    def _normalize_angles(self):
        """ DEPRECATED FUNCTION??? SEEMS TO NOT BE USED ANYWHERE
        """
        self.horizontal_angle = self.horizontal_angle % 360.0
        if self.horizontal_angle<0:
            self.horizontal_angle += 360.0
        if self.vertical_angle>self.max_vertical_angle:
            self.vertical_angle = self.max_vertical_angle
        elif self.vertical_angle<-self.max_vertical_angle:
            self.vertical_angle = -self.max_vertical_angle
    
    def add_orientation_angles(self, h_angle, v_angle):
        """ DEPRECATED FUNCTION??? SEEMS TO NOT BE USED ANYWHERE
        """
        self.horizontal_angle += h_angle
        self.vertical_angle += v_angle
        self._normalize_angles()
        return True
    
    def look_at(self, target):
        """ DEPRECATED FUNCTION??? SEEMS TO NOT BE USED ANYWHERE
        """
        position = self.get_position()
        assert(position[0]!=target[0] and
               position[1]!=target[1] and
               position[2]!=target[2])
        direction = target - position
        direction /= np.linalg.norm(direction)
        self.vertical_angle = -np.arcsin(direction[1])*180/np.pi
        self.horizontal_angle = -(np.arctan2(-direction[0], -direction[2])*180/np.pi)
        self._normalize_angles()
        return True
    
    def get_proj_view_matrix(self):
        """ Returns:
                proview -- a numpy array of size 16 and ndim of 2, containing
                           the information for the camera's projection-view
                           matrix.
        """
        proview = mop.my_glMultiplyMatricesf(self.get_projection_matrix(), self.get_view_matrix())
        return proview
    
    def set_view_matrix(self, new_view_matrix):
        """ Sets a new matrix as view matrix for the camera.
            
            Input parameters:
                new_view_matrix -- a numpy array of size 16 and ndim of 2
                                   containing the new view matrix.
            
            Returns:
                True
        """
        self.view_matrix = new_view_matrix
        return True
    
    def set_projection_matrix(self, new_proj_matrix):
        """ Sets a new matrix as projection matrix for the camera.
            
            Input parameters:
                new_proj_matrix -- a numpy array of size 16 and ndim of 2
                                   containing the new projection matrix.
            
            Returns:
                True
        """
        self.projection_matrix = new_proj_matrix
        return True
    
    def update_projection(self):
        """ Updates the projection matrix. If you change any parameter of the
            camera, you should use this function right after so the changes can
            be applied.
            
            Returns:
                True
        """
        self.projection_matrix = self._get_projection_matrix()
        return True
    
    def update_fog(self):
        """ Updates automatically the fog. This function was created to avoid
            errors when the fog values were changed manually, this way, the fog
            start and end will always have constant values. If you want to
            expand or decrease the fog distance, change the self.min_zfar
            instead of changing this function.
            
            Returns:
                True
        """
        self.fog_end = self.z_far
        self.fog_start = self.fog_end - self.min_zfar
        return True
    
    def print_params(self):
        """ Prints camera parameters in the terminal. Method created only for
            debugging purposes. It will come out in the final distribution?
        """
        logger.debug("######## GLCAMERA PARAMETERS ########")
        logger.debug("<= z_near    => {}".format(self.z_near))
        logger.debug("<= z_far     => {}".format(self.z_far))
        logger.debug("<= fog_start => {}".format(self.fog_start))
        logger.debug("<= fog_end   => {}".format(self.fog_end))
        logger.debug("<= position  => {}".format(self.get_position()))
        logger.debug("######## GLCAMERA PARAMETERS ########")
        return True
    
    def print_matrices(self):
        """ Prints camera matrices in the terminal. Method created only for
            debugging purposes. It will come out in the final distribution?
        """
        logger.debug("######## GLCAMERA MATRICES ########")
        logger.debug("<= view_matrix => {}".format(self.view_matrix))
        logger.debug("<= projection_matrix => {}".format(self.projection_matrix))
        logger.debug("######## GLCAMERA MATRICES ########")
        return True
    
