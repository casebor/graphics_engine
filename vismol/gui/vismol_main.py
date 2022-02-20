#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#

import gi, sys
gi.require_version("Gtk", "3.0")
from gi.repository import Gtk, Gdk
import os
import GTKGUI
gladepath = os.path.split(GTKGUI.__file__)[:-1]
gladepath = os.path.join(*gladepath, "MainWindow.glade")
from GTKGUI.gtkWidgets.filechooser import FileChooser

class VismolMainWindow():
    """ Class doc """
    
    def __init__(self, vismol_session=None, filein=None):
        """ Class initialiser """
        self.builder = Gtk.Builder()
        self.builder.add_from_file(gladepath)
        self.builder.connect_signals(self)
        self.window = self.builder.get_object("MainWindow")
        self.window.set_default_size(800, 600)
        # Status Bar
        self.statusbar_main = self.builder.get_object("MainStatusBar")
        self.statusbar_main.push(1,"Welcome to VisMol")
        self.paned_V = self.builder.get_object("paned_V")
        self.vismol_session = vismol_session
        self.vismol_session.main_session = self
        self.window.connect("key-press-event", self.vismol_session.glwidget.key_pressed)
        self.window.connect("key-release-event", self.vismol_session.glwidget.key_released)
        self.menu_box = self.builder.get_object("toolbutton16")
        self.box2 = self.builder.get_object("box2")
        self.selection_box = self.vismol_session.selection_box
        self.menu_box.add(self.selection_box)
        #-------------------------------------------------------------------      
        #                         notebook_V1
        #-------------------------------------------------------------------
        #self.notebook_V1 = Gtk.Notebook()
        #print (self.notebook_V1.set_tab_pos(Gtk.PositionType.LEFT))
        #self.page1 = Gtk.Box()
        #self.page1.set_border_width(5)
        
        #self.text_view = Gtk.TextView()
        #self.text_view.set_editable(True)
        #self.page1.add( self.text_view)
        
        #self.page1.add(Gtk.Label("Here is the content of the first section."))
        #self.notebook_V1.append_page(self.page1, Gtk.Label("Logs"))
        
        #-------------------------------------------------------------------      
        #                         notebook_H1
        #-------------------------------------------------------------------
        self.notebook_H1 = Gtk.Notebook()
        self.ScrolledWindow_notebook_H1 = Gtk.ScrolledWindow()
        self.treeview = GtkMainTreeView(vismol_session)
        self.ScrolledWindow_notebook_H1.add(self.treeview)
        self.notebook_H1.append_page(self.ScrolledWindow_notebook_H1, Gtk.Label("Objects"))
        # the label we use to show the selection
        self.label = Gtk.Label()
        self.label.set_text("")
        
        #-------------------------------------------------------------------
        #                         notebook_H2
        #-------------------------------------------------------------------
        self.notebook_H2 = Gtk.Notebook()
        #-------------------------------------------------------------------
        self.paned_H = Gtk.Paned(orientation = Gtk.Orientation.HORIZONTAL)
        self.button = Gtk.Button(label="Click Here")
        #-------------------------------------------------------------------
        self.vismol_session = vismol_session
        self.filechooser   = FileChooser()
        #-------------------------------------------------------------------
        
        self.container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        self.command_line_entry = Gtk.Entry()
        
        if self.vismol_session is not None:
            #player
            self.container.pack_start(self.vismol_session.glwidget, True, True, 0)
            self.traj_frame = self.vismol_session.trajectory_frame
            self.notebook_H2.append_page(self.container, Gtk.Label("view"))
            self.notebook_H2.append_page(Gtk.TextView(), Gtk.Label("logs"))
            self.HBOX = Gtk.Box(orientation = Gtk.Orientation.VERTICAL, spacing = 0)
            self.HBOX.pack_start(self.notebook_H1, True, True, 0)
            self.HBOX.pack_start(self.traj_frame, False, False, 1)
            self.paned_H.add(self.HBOX)
            self.paned_H.add(self.notebook_H2)
            self.paned_V.add(self.paned_H)
        self.window.connect("delete-event",    Gtk.main_quit)
        self.window.show_all()
        
        if filein:
            self.vismol_session.load(filein)
        Gtk.main()
    
    def gtk_load_files(self, button):
        filename = self.filechooser.open()
        if filename:
            self.vismol_session.load(filename)
    
    def run(self):
        """ Function doc """
        Gtk.main()
    
    def menubar_togglebutton1(self, button):
        """ Function doc """
        if button.get_active():
            state = "on"
            self.vismol_session._picking_selection_mode = True
            button.set_label("Picking")
        else:
            state = "off"
            self.vismol_session._picking_selection_mode = False
            button.set_label("Viewing")
        print("was turned", state)
    
    def test (self, widget):
        """ Function doc """
        container = Gtk.Box (orientation = Gtk.Orientation.VERTICAL)
        container.pack_start(self.notebook_V1, True, True, 0)
        container.pack_start(self.command_line_entry, False, False, 0)
        self.paned_V.add(container)
        self.paned_V.show()
        self.window.show_all()
    
    def on_toolbutton_trajectory_tool (self, button):
        """ Function doc """
        print (button)


class GtkMainTreeView(Gtk.TreeView):
    """ Class doc """
    
    def __init__ (self, vismol_session):
        """ Class initialiser """
        Gtk.TreeView.__init__(self)
        self.vismol_session = vismol_session
        self.treeview_menu = TreeViewMenu(self)
        self.store = vismol_session.Vismol_Objects_ListStore
        self.set_model(self.store)
        #----------------------------------------------------------------------
        # the cellrenderer for the second column - boolean rendered as a toggle
        renderer_toggle = Gtk.CellRendererToggle()
        # the second column is created
        column_in_out = Gtk.TreeViewColumn("", renderer_toggle, active=0)
        # and it is appended to the treeview
        self.append_column(column_in_out)
        # connect the cellrenderertoggle with a callback function
        renderer_toggle.connect("toggled", self.on_toggled)
        # the cellrenderer for text columns
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("id", renderer_text, text=1)
        self.append_column(column)
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("Object", renderer_text, text=2)
        self.append_column(column)
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("Atoms",  renderer_text, text=3)
        self.append_column(column)
        renderer_text = Gtk.CellRendererText()
        column = Gtk.TreeViewColumn("Frames", renderer_text, text=4)
        self.append_column(column)
        self.connect("button-release-event", self.on_treeview_Objects_button_release_event )
        #----------------------------------------------------------------------
    
    def on_toggled(self, widget, path):
        # the boolean value of the selected row
        current_value = self.store[path][0]
        # change the boolean value of the selected row in the model
        self.store[path][0] = not current_value
        if self.store[path][0]:
            obj_index = self.store[path][1]
            self.vismol_session.enable_by_index(int(obj_index))
            self.vismol_session.glwidget.queue_draw()
        else:
            obj_index = self.store[path][1]
            self.vismol_session.disable_by_index(int(obj_index))
            self.vismol_session.glwidget.queue_draw()
    
    def append(self, visObj):
        i = visObj.index
        data = [visObj.active         , 
               str(i)                 ,
               visObj.name            , 
               str(len(visObj.atoms)) , 
               str(len(visObj.frames)),
               ]
        print (data)
        self.store.append(data)
    
    def remove (self):
        """ Function doc """
        pass
        
    def on_treeview_Objects_button_release_event(self, tree, event):
        if event.button == 3:
            selection = self.get_selection()
            model = self.get_model()
            (model, iter) = selection.get_selected()
            if iter != None:
                self.selectedID  = str(model.get_value(iter, 1))
                self.selectedObj = str(model.get_value(iter, 2))
                self.treeview_menu.open_menu(self.selectedObj)
        
        if event.button == 2:
            selection= tree.get_selection()
            model= tree.get_model()
            (model, iter) = selection.get_selected()
            self.selectedID  = int(model.get_value(iter, 1))
            visObj = self.vismol_session.vismol_objects_dic[self.selectedID]
            self.vismol_session.center(visObj)
        
        if event.button == 1:
            print ("event.button == 1:")


class TreeViewMenu:
    """ Class doc """
    
    def __init__ (self, treeview):
        """ Class initialiser """
        pass
        self.treeview = treeview
        self.filechooser   = FileChooser()
        functions = {
                "test":self.f1 ,
                "f1": self.f1,
                "f2": self.f2,
                "delete": self.f3,
        }
        self.build_glmenu(functions)
    
    def f1 (self, visObj = None ):
        """ Function doc """
        selection        = self.treeview.get_selection()
        (model, iter)    = selection.get_selected()
        self.selectedID  = int(model.get_value(iter, 1))
        visObj = self.treeview.vismol_session.vismol_objects_dic[self.selectedID]
        infile = self.filechooser.open()
        self.treeview.vismol_session.load_xyz_coords_to_vismol_object(infile , visObj)
        model[iter][4] = str(len(visObj.frames))
        print("f1")
    
    def f2 (self, visObj = None):
        """ Function doc """
        self.treeview.vismol_session.go_to_atom_window.OpenWindow()
        print("f2")
    
    def f3 (self, visObj = None):
        """ Function doc """
        selection     = self.treeview.get_selection()
        (model, iter) = selection.get_selected()
        self.selectedID  = int(model.get_value(iter, 1))
        del self.treeview.vismol_session.vismol_objects_dic[self.selectedID]
        self.treeview.store .clear()
        for vobj_index, vis_object in self.treeview.vismol_session.vismol_objects_dic.items():
            print ("\n\n",vis_object.name,"\n\n")
            data = [vis_object.active          , 
                    str(vis_object.index),
                    vis_object.name            , 
                    str(len(vis_object.atoms)) , 
                    str(len(vis_object.frames)),
                   ]
            model.append(data)
        self.treeview.vismol_session.glwidget.queue_draw()
        print("f3")
    
    def build_glmenu (self, menu_items = None):
        """ Function doc """
        self.glMenu = Gtk.Menu()
        for label in menu_items:
            mitem = Gtk.MenuItem(label)
            mitem.connect("activate", menu_items[label])
            self.glMenu.append(mitem)
            mitem = Gtk.SeparatorMenuItem()
            self.glMenu.append(mitem)
        self.glMenu.show_all()
    
    def open_menu (self, visObj = None):
        """ Function doc """
        print (visObj)
        self.glMenu.popup(None, None, None, None, 0, 0)


class GtkMainTreeView_old(Gtk.TreeView):
    """ Class doc """
    
    def __init__(self, vismol_session):
        """ Function doc """
        Gtk.TreeView.__init__(self)
        self.connect("button-release-event", self.on_treeview_Objects_button_release_event)
        self.vismol_session = vismol_session
        columns = [" " ,
                   "id",
                   "Object",
                   "Atoms" ,
                   "Frames"]
        self.liststore  = Gtk.ListStore(bool,str ,str, str, str)
        for i, column in enumerate(columns):
            if i == 0:
                cell = Gtk.CellRendererToggle()
                cell.set_property("activatable", True)
                cell.connect("toggled", self.on_chk_renderer_toggled, self.liststore)
                col = Gtk.TreeViewColumn(None, cell )
                col.add_attribute(cell, "active", 0)
            else:
                cell = Gtk.CellRendererText()
                col = Gtk.TreeViewColumn(column, cell, text=i)
            self.append_column(col)
        self.treeview_menu = TreeViewMenu()
    
    def on_chk_renderer_toggled(self, cell, path, model):
        """ Function doc """
        model[path][0] = not model[path][0]
        print (model[path][0], model, )
        for item in model:
            print (item)
    
    def append(self, visObj):
        """ Function doc """
        i = self.vismol_session.vismol_objects.index(visObj)
        data = [visObj.active, str(i), visObj.name, str(len(visObj.atoms)), str(len(visObj.frames))]
        self.liststore.append(data)
        self.set_model(self.liststore)
    
    def refresh_gtk_main_treeview (self):
        """ Function doc """
        model = self.liststore  
        model.clear()
        n = 0
        i = 1
        for vis_object in self.vismol_session.vismol_objects:
            data = [vis_object.active, str(i), vis_object.name,
                    str(len(vis_object.atoms)), str(len(vis_object.frames))]
            model.append(data)
            i +=1
            n = n + 1
        self.set_model(model)
    
    def on_treeview_Objects_button_release_event(self, tree, event):
        """ Function doc """
        if event.button == 3:
            selection = self.get_selection()
            model = self.get_model()
            (model, iter) = selection.get_selected()
            if iter != None:
                self.selectedID  = str(model.get_value(iter, 1))
                self.selectedObj = str(model.get_value(iter, 2))
                self.treeview_menu.glMenu.popup(None, None, None, None, 0, 0)
        
        if event.button == 2:
            self.refresh_gtk_main_self.treeView()
            print ("button == 2")
        
        if event.button == 1:
            print ("event.button == 1:")
    
    def on_treemenu_item_selection (self, widget, event=None, data=None):
        """ Function doc """
        if widget == self.builder.get_object("menuitem5_rename"):
            tree = self.builder.get_object("treeview1")
            selection = tree.get_selection()
            model = tree.get_model()
            (model, iter) = selection.get_selected()
            obj_index = model.get_value(iter, 1)
            self.vismol_session.edit_by_index(int(obj_index)-1)
            self.vismol_session.glwidget.vm_widget.editing_mols = not self.vismol_session.glwidget.vm_widget.editing_mols
        tree = self.builder.get_object("treeview1")
        selection = tree.get_selection()
        model = tree.get_model()
        (model, iter) = selection.get_selected()
        obj_index = model.get_value(iter, 1)
        visObj = self.vismol_session.vismol_objects[(int(obj_index)-1)]
        
        if widget == self.builder.get_object("menuitem_center"):
            self.vismol_session.glwidget.vm_widget.center_on_coordinates(visObj, visObj.mass_center)
        
        if widget == self.builder.get_object("menu_show_lines"):
            visObj.lines_actived = True
        
        if widget == self.builder.get_object("menu_show_sticks"):
            visObj.sticks_actived = True
        
        if widget == self.builder.get_object("menu_show_spheres"):
            visObj.spheres_actived = True
        
        if widget == self.builder.get_object("menu_show_ribbons"):
            visObj.ribbons_actived = True
        
        if widget == self.builder.get_object("menu_show_dots"):
            visObj.dots_actived = True
            self.vismol_session.glwidget.vm_widget.queue_draw()
        
        if widget == self.builder.get_object("menu_hide_lines"):
            visObj.lines_actived = False
        
        if widget == self.builder.get_object("menu_hide_sticks"):
            visObj.sticks_actived = False
        
        if widget == self.builder.get_object("menu_hide_spheres"):
            visObj.spheres_actived = False
        
        if widget == self.builder.get_object("menu_hide_ribbons"):
            visObj.ribbons_actived = False
            
        if widget == self.builder.get_object("menu_hide_dots"):
            visObj.dots_actived = False
            self.vismol_session.glwidget.vm_widget.queue_draw()
