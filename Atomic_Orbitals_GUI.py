# Developed by Dr. Matthew D. Hanson 
# This module can be used to plot 3-D representations of the orbitals used
# by chemists to describe the electronic states of the hydrogen atom or
# hydrogenlike ions. The user need only specify the n,l, and m quantum
# numbers of the desired orbitals, where a choice of m greater than 0 will
# result in two orbitals from the real and imaginary parts of psi_nlm. Note
# that radial wavefunctions are provided to aid with interpretation. Also 
# note that intrinsic limitations of plotting 3-D graphics may cause some 
# artifacts to appear when rotating the orbitals. The portions of this code
# are sectioned off in the code blocks that appear in the Jupyter notebook. 
# Block 1 #
import math
import numpy                #requires separate installation
import scipy                #requires separate installation
from scipy.optimize import root
from scipy import special
import sympy as sym         #requires separate installation
import tkinter as tk        #requires separate installation
import ast
import matplotlib as mpl    #requires separate installation
import matplotlib.pyplot as plt 
from matplotlib.pyplot import cm
from tvtk.util import ctf   #requires separate installation
from tvtk.util.ctf import ColorTransferFunction     #used for hacky colormap editing
from tvtk.util.ctf import PiecewiseFunction         #used for hacky opacity editing
print('Initial Packages successfully imported.')    #prints out a message telling the user this section ran successfully

# Block 2 #
def submit_button():  #submits updated values for all the variables
    global entries,checks,drops
    print('Variable values set:')
    entries = [entry_variables[i].get() for i in range(len(entry_variables))]
    checks = [check_variables[i].get() for i in range(len(check_variables))]
    drops = [drop_variables[i].get() for i in range(len(drop_variables))]
    print('You may close the window now.')

def defaults_button():  #resets all the values to the defaults
    for i in range(len(entry_variables)):
        entry_variables[i].set(default_entry_variables[i])
    for i in range(len(check_variables)):
        check_variables[i].set(default_check_variables[i])
    for i in range(len(drop_variables)):
        drop_variables[i].set(default_drop_variables[i])
    print('All variables reset to default values:')

# Data boxes to hold information about the GUI; ordered label, box type, row, column
Data_boxes = numpy.array([['Principal Quantum Number n: ','entry',0,0], ['Allowed Values: 1,2,3,...','label',0,2],
                          ['Angular Momentum Quantum Number l: ','entry',1,0], ['Allowed Values: 0,1,...,n-1','label',1,2],
                          ['Magnetic Quantum Number m: ','entry',2,0], ['Allowed Values: -l,-l+1,...,l-1,l','label',2,2],
                          ['Colormap: ','drop',3,0],['Reverse Colormap: ','check',3,3],
                          ['Colormap points: ','entry',4,0], ['Allowed Values: 1,3,5,...,255','label',4,2], ['Cut Quadrant for m = 0 Orbitals: ','check',4,3],
                          ['Plot Mode: ','drop',5,0],['Scale Wavefunction By Maximum Value: ','check',5,3],
                          ['Fractional Cutoff Value: ','entry',7,0], ['Allowed Values in [0,1]','label',7,2],
                          ['Plot All Nodes: ','check',8,0],['Nodal Surface Opacity: ','entry',8,3], ['Allowed Values in [0,1]','label',8,5],
                          ['Plot Radial Nodes: ','check',9,0],['Radial Node RGB Values: ','entry',9,3], ['Allowed Values in [0,1]','label',9,5],
                          ['Plot \u03b8 Nodes: ','check',10,0],['\u03b8 Node RGB Values: ','entry',10,3], ['Allowed Values in [0,1]','label',10,5],
                          ['Plot \u03c6 Nodes: ','check',11,0],['\u03c6 Node RGB Values: ','entry',11,3], ['Allowed Values in [0,1]','label',11,5],
                          ['RGBs Are Tuples','label',12,4],
                          ['Opacity Scaling Exponent: ','entry',13,0],['Must Be 0 or Above','label',13,2],['Maximum Opacity: ','entry',13,3], ['Allowed Values in [0,1]','label',13,5],
                          ['Opacity Shift: ','entry',14,0],['Allowed Values in [0,1]','label',14,2],
                          ['Use White Background: ','check',15,0],['Outline Box: ','check',15,3],
                          ['Round Axes to Multiples of 10 a.u.: ','check',16,0],
                          ['Maximum Radial Distance (a.u.) Per n: ','drop',17,0],['Number of Grid Points Per n: ','drop',17,3],
                          ['Defaults','defaults_button',18,2],['Submit','submit_button',18,3]])
# Options ordered in the same order they appear in Data_boxes; not all available colormaps are actually given #
Drop_options = [['seismic','coolwarm','bwr','jet','viridis','plasma','cool','hot','rainbow'],
                ['orbital','probability density','original'],
                [str(2*i+6) for i in range(13)],
                [str(5*i+10) for i in range(10)]]
default_entry_variables = ['2','1','0','21','0.0','1.0','(0,1,0)','(0,1,0)','(0,1,0)','1.0','0.6','0.0']
default_check_variables = [False,False,True,False,False,False,False,True,False,False]
default_drop_variables = ['seismic','orbital','12','20']

All_valid = False  #master switch that needs to be True to proceed
while All_valid == False:  #the user can't escape the GUI sequence until all the variables are valid
   # Create a new window with the title 'Atomic Orbital Setup'
    window = tk.Tk()
    window.title("Atomic Orbital Setup")
    frame = tk.Frame(relief=tk.SUNKEN, borderwidth=3)
    frame.pack()

    # Setting up the GUI variables and result arrays with their defaults # 
    entry_variables = [tk.StringVar() for i in range(len(default_entry_variables))]
    check_variables = [tk.BooleanVar() for i in range(len(default_check_variables))]
    drop_variables = [tk.StringVar() for i in range(len(default_drop_variables))]
    drop_widths = [5,5,5,5,5]  #how wide each drop box will be
    for i in range(len(entry_variables)):
        entry_variables[i].set(default_entry_variables[i])
    for i in range(len(check_variables)):
        check_variables[i].set(default_check_variables[i])
    for i in range(len(drop_variables)):
        drop_variables[i].set(default_drop_variables[i])
    entries = [entry_variables[i].get() for i in range(len(entry_variables))]
    checks = [check_variables[i].get() for i in range(len(check_variables))]
    drops = [drop_variables[i].get() for i in range(len(drop_variables))]
    
    # Create the widgets for the GUI settings #
    entry_counter = -1  #counters to help us sort our variables and widgets
    check_counter = -1
    drop_counter = -1
    for idx, data in enumerate(Data_boxes):  #makes labeled widgets of various types based on Data_boxes
        if data[1] == 'entry':  #sets up the entry boxes
            entry_counter += 1
            label = tk.Label(master=frame, text=data[0])
            obj = tk.Entry(master=frame, width=10, textvariable = entry_variables[entry_counter])
            label.grid(row=data[2], column=data[3], sticky="e")
            obj.grid(row=int(data[2]), column=int(data[3])+1,padx=5,pady=2)
        if data[1] == 'check':  #sets up the check boxes
            check_counter += 1
            label = tk.Label(master=frame, text=data[0])
            obj = tk.Checkbutton(master=frame, width=5, variable = check_variables[check_counter],onvalue=True, offvalue=False)
            label.grid(row=data[2], column=data[3], sticky="e",padx=5,pady=2)
            obj.grid(row=int(data[2]), column=int(data[3])+1)
        if data[1] == 'drop':  #sets up the drop down menus
            drop_counter += 1
            label = tk.Label(master=frame, text=data[0])
            obj = tk.OptionMenu(frame, drop_variables[drop_counter], *Drop_options[drop_counter])
            label.grid(row=data[2], column=data[3], sticky="e")
            obj.grid(row=int(data[2]), column=int(data[3])+1,padx=5,pady=2)
            obj.config(fg='black', bg='white', width = drop_widths[drop_counter])
        if data[1] == 'label':  #sets up labels with no other widgets
            label = tk.Label(master=frame, text=data[0])
            label.grid(row=data[2], column=data[3], sticky="w",padx=5,pady=2)
        if data[1] == 'submit_button':  #button to submit the data
            obj = tk.Button(master=frame, text=data[0], fg='black', bg='white', command=submit_button)
            obj.grid(row=data[2], column=data[3], ipadx=10, sticky='w')
        if data[1] == 'defaults_button':  #button to reset things to defaults
            obj = tk.Button(master=frame, text=data[0], fg='black', bg='white', command=defaults_button)
            obj.grid(row=data[2], column=data[3], ipadx=10, sticky='e')

    window.mainloop()  #activates the GUI

    # Validating the user-defined data to ensure the program will not break #
    QN_valid = False  #tells if quantum numbers are valid
    try:  #Testing the quantum numbers
        n = int(entries[0])  #principal quantum number
        l = int(entries[1])  #angular momentum quantum number
        m = int(entries[2])  #magnetic quantum number
        if n > 0 and l >= 0 and l <= n-1 and m >= -l and m <= l:
            QN_valid = True
        else:  #Tells how to actually get valid quantum numbers
            print('The quantum numbers must be integers that obey the appropriate ranges.')
    except:
        print('The quantum numbers must be integers that obey the appropriate ranges.')
    color_points_valid = False  #tells if number of colormap points is valid
    try:
        colormap_points = int(entries[3])   #number of points used to change colors; more gives smoother colors, max 255
        if colormap_points > 0:  #we will fix even values later
            color_points_valid = True
        else:
            print('Please provide a valid number of points.')
    except:
        print('Please provide a valid number of points.')
    cutoff_valid = False
    try:    
        cutoff = float(entries[4])   #wavefunction magnitude fractional cutoff for plotting
        if cutoff >= 0.0 and cutoff <= 1.0:
            cutoff_valid = True
        else:
            print('Cutoff must be in range [0,1].')
    except:
        print('Cutoff must be in range [0,1].')
    node_opacity_valid = False
    try:    
        node_opacity = float(entries[5])   #opacity of the nodal surfaces
        if node_opacity >= 0.0 and node_opacity <= 1.0:
            node_opacity_valid = True
        else:
            print('Node opacity must be in range [0,1].')
    except:
        print('Node opacity must be in range [0,1].')
    RGB_valid = False
    length_valid = False
    try:
        radial_node_color = tuple(ast.literal_eval(entries[6]))  #RGB tuples
        theta_node_color = tuple(ast.literal_eval(entries[7]))
        phi_node_color = tuple(ast.literal_eval(entries[8]))
        if len(radial_node_color) == 3 and len(theta_node_color) == 3 and len(theta_node_color) == 3:
            length_valid = True
        rad_test = [(radial_node_color[i] >= 0 and radial_node_color[i] <= 1) for i in range(3)]
        theta_test = [(theta_node_color[i] >= 0 and theta_node_color[i] <= 1) for i in range(3)]
        phi_test = [(phi_node_color[i] >= 0 and phi_node_color[i] <= 1) for i in range(3)]
        if rad_test and theta_test and phi_test and length_valid:
            RGB_valid = True
        else:
            print('RGB values must be length 3 tuples with values between 0 and 1')
    except:
        print('RGB values must be length 3 tuples with values between 0 and 1')
    opacity_exp_valid = False
    try:
        opacity_exponent = float(entries[9])  #used to scale the opacity; higher value makes low amplitude more transparent
        if opacity_exponent > 0:
            opacity_exp_valid = True
        else:
            print('Opacity exponent should be greater than 0.')
    except:
        print('Opacity exponent should be greater than 0.')
    opacity_max_valid = False
    try:
        opacity_factor = float(entries[10])  #maximum value of the opacity
        if opacity_factor >= 0 and opacity_factor <= 1:
            opacity_max_valid = True
        else:
            print('Opacity limit should be in the range [0,1].')
    except:
        print('Opacity limit should be in the range [0,1].')
    opacity_shift_valid = False
    try:
        opacity_shift = float(entries[11])  #gives the minimum of the opacity. you should use a non-0 cutoff if this is not set to 0
        if opacity_shift >= 0 and opacity_shift <= 1:
            opacity_shift_valid = True
        else:
            print('Opacity shift should be in the range [0,1].')
    except:
        print('Opacity shift should be in the range [0,1].')

    reverse_colormap = checks[0]  #True reverses colormap, False leaves it alone
    m0_cut = checks[1]      #True will make orbitals with m=0 orbitals have the first quadrant cut out for viewing, False leaves it alone
    if checks[2]:
        wavefunction_mode = 'scaled'  #scaled divides wavefunction by its maximum value, actual gives the unscaled values
    if checks[2] == False:
        wavefunction_mode = 'actual'
    plot_nodes = checks[3]          #if True, all nodal surfaces will be plotted for the orbitals, if False checks will be made for each one
    plot_radial_nodes = checks[4]   #if True, radial nodal surfaces will be plotted, if False, they may not. Overrules plot_nodes = False
    plot_theta_nodes = checks[5]    #if True, theta nodal surfaces will be plotted, if False, they may not. Overrules plot_nodes = False
    plot_phi_nodes = checks[6]      #if True, phi nodal surfaces will be plotted, if False, they may not. Overrules plot_nodes = False
    if checks[7]:
        light_mode = 'white'    #used to determine if the plot background is white or black
    if checks[7] == False:
        light_mode = 'black'
    box_outline = checks[8]     #if True plots a box around the edges of the grid, if False does not
    nice_labels = checks[9]    #if True makes axis labels a multiple of 10, if False gives real values

    colormap = drops[0]
    plot_mode = drops[1] 
    Node_mode = 'mesh'
    rmax = int(drops[2])*n    #Maximum distance along x,y,z (units of a_0 or a.u.). The default is sufficient for most applications
    N = int(drops[3])*n       #number of grid points along r,x,y,z axes used throughout the program (~35*N gives high quality)

    if QN_valid and color_points_valid and cutoff_valid and node_opacity_valid and RGB_valid and opacity_exp_valid and opacity_max_valid and opacity_shift_valid:
        All_valid = True  #terminates the while loop and allows the program to continue
        print('Variables successfully defined.')
    else:
        print('Not all chosen values are valid. Try again.')

# Block 3 #
from mayavi import mlab     #requires separate installation
print('mlab imported from mayavi.')

# Block 4 #
# Console output detailing the settings #
if m == 0 or (m != 0 and plot_mode == 'original'):
    String = 'Orbital requested for n=' + str(n) + ', l=' + str(l) + ', m=' + str(abs(m))
if m != 0 and plot_mode != 'original':
    String = 'Orbital requested for n=' + str(n) + ', l=' + str(l) + ', m=' + str(abs(m)) + ',-' + str(abs(m))
orbital_names = ['s','p','d','f','g','h']
if l < 6:
    if reverse_colormap:
        String += ' ('  + str(n) + orbital_names[l] + ') will be plotted in ' + plot_mode + ' mode with the ' + colormap + '_r colormap.'
    else:
        String += ' ('  + str(n) + orbital_names[l] + ') will be plotted in ' + plot_mode + ' mode with the ' + colormap + ' colormap.'
if l >= 6:
    if reverse_colormap:
        String += ' will be plotted in ' + plot_mode + ' mode with the ' + colormap + '_r colormap.'
    else:
        String += ' will be plotted in ' + plot_mode + ' mode with the ' + colormap + ' colormap.'
print(String)
if colormap_points%2 == 1:
    print('Volumetric colormap will use ' + str(colormap_points) + ' points.')
if colormap_points%2 == 0:
    print('An even number of colormaps points was supplied: ' + str(colormap_points) + '. An odd number of points will be used instead: ' + str(colormap_points+1))
    colormap_points += 1
if plot_nodes or (plot_radial_nodes and plot_theta_nodes and plot_phi_nodes):
    print('Nodes will be plotted for the radial and angular wavefunctions with an opacity of ' + str(node_opacity) + '.')
if plot_nodes == False and (plot_radial_nodes or plot_theta_nodes or plot_phi_nodes) and (plot_radial_nodes and plot_theta_nodes and plot_phi_nodes) == False:
    print('Nodes will be plotted for portions of the wavefunction with a mesh opacity of ' + str(node_opacity) + '.')
if wavefunction_mode == 'scaled':
    if (m0_cut and m == 0):
        print('The orbital wavefunction will be scaled by its maximum value and the first quadrant will be cut out.')
    if m0_cut == False or m != 0:
        print('The orbital wavefunction will be scaled by its maximum value.')
if wavefunction_mode == 'actual':
    if (m0_cut and m == 0):
        print('The orbital wavefunction will be unmodified and the first quadrant will be cut out.')
    if m0_cut == False or m != 0:
        print('The orbital wavefunction will be unmodified.')
print('Volumetric data will be displayed with a maximum opacity of ' + str(opacity_factor) + ' and opacity scaling power of ' + str(opacity_exponent) + '.')
print('The grid will span -' + str(rmax) + ' a.u. to ' + str(rmax) + ' a.u. along x,y, and z and will use ' + str(N) + ' grid points along x,y, and z.')
print('This will generate a total of ' + str(N**3) + ' grid points for the wavefunction.')
if N**3 > 10**7 and N**3 < 10**8:
    print('### Note that this is a large grid and some waiting is to be expected. ###')
if N**3 > 10**8:
    print('### Note that this grid is extremely large and the program is unlikely to finish running. ###')

# Block 5 #
# Defining a class and function for setting colormaps and opacity #
class MplColorHelper:   #used to get the RGB information from a desired colormap
    def __init__(self, cmap_name, start_val, stop_val):
        self.cmap_name = cmap_name
        self.cmap = plt.get_cmap(cmap_name)
        self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)
    def get_rgb(self, val):   #used to turn the data values along a colormap into RGB values
        return self.scalarMap.to_rgba(val)

def Modify_Color_Opacity(volume,orbital,magnitude): 
    # Modifying the orbital colormap #
    orbital_max = numpy.amax(numpy.abs(orbital))
    c = ColorTransferFunction()  #sets up a blank ctf 
    if magnitude == False:
        values = numpy.linspace(-orbital_max,orbital_max,colormap_points)   #gets the values for the custom colormap
        for i in range(colormap_points):
            c.add_rgb_point(values[i], color_array[i][0], color_array[i][1], color_array[i][2]) #Add enough points to CTF to make a smooth colormap
    if magnitude:  #only uses half the colormap for the magnitude plot
        values = numpy.linspace(0,orbital_max,colormap_points//2+1)   #gets the values for the custom colormap
        for i in range(colormap_points//2+1):
            c.add_rgb_point(values[i], color_array[colormap_points//2+i][0], color_array[colormap_points//2+i][1], color_array[colormap_points//2+i][2]) 
    volume._volume_property.set_color(c)    #sets the new colormap
    
    # Modifying the orbital opacity #
    otf = PiecewiseFunction()               #sets up a blank opacity transfer function
    values = numpy.linspace(-orbital_max,orbital_max,colormap_points)  #set up like the color array, but only positive values
    for i in range(colormap_points):    #scales the opacity of points exponentially for better visuals of interior structure
        otf.add_point(values[i],opacity_factor*((numpy.abs(-colormap_points//2+i+1)/(colormap_points//2))**opacity_exponent))
    volume.volume_property.set_scalar_opacity(otf)  #setting up the new opacity function
    
    # Hacking the colorbar by adding a transparent small surface #
    xhack,yhack,zhack = numpy.mgrid[-1:1:2j, -1:1:2j, -1:1:2j]
    Phihack = numpy.zeros(numpy.shape(xhack))
    if wavefunction_mode == 'scaled':
        Phihack[0][0][0] = -1
        Phihack[1][1][1] = 1
    if wavefunction_mode == 'actual':
        Phihack[0][0][0] = -orbital_max
        Phihack[1][1][1] = orbital_max
    hack = mlab.points3d(xhack, yhack, zhack, Phihack, opacity=0, colormap=colormap)    #sets up a completely transparent object with the right colormap
    if reverse_colormap == True:
        hack.module_manager.scalar_lut_manager.reverse_lut = True
    
    #setting the colorbar with the transparent object
    if wavefunction_mode == 'scaled' and (plot_mode == 'original' or plot_mode == 'orbital') and magnitude == False:
        colorbar = mlab.colorbar(object = hack, title='Wavefunction Value (Scaled)', orientation='horizontal', nb_labels=5, label_fmt='%.3f') 
    if wavefunction_mode == 'scaled' and plot_mode == 'probability density' and magnitude == False:
        colorbar = mlab.colorbar(object = hack, title='Probability Density (Scaled)', orientation='horizontal', nb_labels=5, label_fmt='%.3f') 
    if wavefunction_mode == 'scaled' and magnitude:
        colorbar = mlab.colorbar(object = hack, title='Orbital Magnitude (Scaled)', orientation='horizontal', nb_labels=5, label_fmt='%.3f') 
    if wavefunction_mode == 'actual' and (plot_mode == 'original' or plot_mode == 'orbital') and magnitude == False:
        colorbar = mlab.colorbar(object = hack, title='Wavefunction Value', orientation='horizontal', nb_labels=5, label_fmt='%.3f') 
    if wavefunction_mode == 'actual' and plot_mode == 'probability density' and magnitude == False:
        colorbar = mlab.colorbar(object = hack, title='Probability Density', orientation='horizontal', nb_labels=5, label_fmt='%.3f') 
    if wavefunction_mode == 'actual' and magnitude:
        colorbar = mlab.colorbar(object = hack, title='Orbital Magnitude (Scaled)', orientation='horizontal', nb_labels=5, label_fmt='%.3f') 
    colorbar.scalar_bar_representation.position = [0.1, 0.85]  #moves the colorbar
    colorbar.scalar_bar.unconstrained_font_size = True    #allows you to actually edit the font size
    colorbar.label_text_property.font_size=20
    colorbar.title_text_property.font_size=20
    volume.update_ctf = True

# Defining a function that generates nodal surfaces for a set of nodal value arrays #
def Mayavi_Nodes(r_nodes, theta_nodes, phi_nodes):
    node_grid = 5 #scales number of grid points used in plotting nodes
    tube_rad = 0.0015*rmax+0.006 #generally makes nice looking node meshes
    if nice_labels:
        distance_cap = rmax-rmax%10  #matches up with the label convention for nice_axes; smaller than rmax
    if nice_labels == False:
        distance_cap = rmax  #just use the real rmax value
    if len(r_nodes) > 0 and (plot_nodes or plot_radial_nodes):  #checks if there are nodes and if the plot is requested 
        for i in range(len(r_nodes)):
            tnodes,pnodes = numpy.mgrid[0:numpy.pi:node_grid*n*1j, 0:2*numpy.pi:2*node_grid*n*1j]  #floor reduces the number of grid points
            xnodes = r_nodes[i]*numpy.sin(tnodes)*numpy.cos(pnodes)
            ynodes = r_nodes[i]*numpy.sin(tnodes)*numpy.sin(pnodes)
            znodes = r_nodes[i]*numpy.cos(tnodes)  #spherical coordiante transformations applied to the grids
            mlab.mesh(xnodes,ynodes,znodes,representation=Node_mode,color = radial_node_color,opacity=node_opacity,tube_radius=tube_rad)
    if len(theta_nodes) > 0 and (plot_nodes or plot_theta_nodes):  #checks if there are nodes and if the plot is requested 
        for i in range(len(theta_nodes)):
            if numpy.round(theta_nodes[i],2) != numpy.round(numpy.pi/2,2):  #non-pi/2 nodes are cones
                tnodes,rnodes = numpy.mgrid[0:2*numpy.pi:node_grid*n*1j, 0:rmax:node_grid*n*1j]  # making cylindrical grid
                xnodes = rnodes*numpy.cos(tnodes)  #cylindrical coordinate transform
                ynodes = rnodes*numpy.sin(tnodes)  #cylindrical coordinate transform
                rnodes[numpy.where(rnodes > distance_cap)] = float('NaN')  #crops values where r is too big for the box
                znodes = 1/math.tan(theta_nodes[i])*rnodes  #cylindrical coordinate transform with factor for theta angle
            if numpy.round(theta_nodes[i],2) == numpy.round(numpy.pi/2,2):  #just in case the solver is somewhat inaccurate
                xnodes,ynodes = numpy.mgrid[-rmax:rmax:(node_grid*n*1j) , -rmax:rmax:(node_grid*n*1j)] #floor gives a smaller grid for the nodes
                znodes = 0*xnodes #for pi/2 theta node you just get the x-y plane
            znodes[numpy.where(numpy.abs(znodes) > distance_cap)] = float('NaN')  #crops values where z is too big for the box
            mlab.mesh(xnodes,ynodes,znodes,representation=Node_mode,color = theta_node_color,opacity=node_opacity,tube_radius=tube_rad)
    if len(phi_nodes) > 0 and (plot_nodes or plot_phi_nodes):  #checks if there are nodes and if the plot is requested 
        for i in range(len(phi_nodes)):
            if numpy.round(phi_nodes[i],2) != 0.0 and numpy.round(phi_nodes[i],2) != numpy.round(numpy.pi/2,2):
                ynodes,znodes = numpy.mgrid[-distance_cap:distance_cap:(node_grid*n*1j) , -distance_cap:distance_cap:(node_grid*n*1j)] #floor gives a smaller grid for the nodes
                xnodes = -math.tan(phi_nodes[i] + numpy.pi/2)*ynodes  #consistent with phi nodes in range (0,pi] only
            if numpy.round(phi_nodes[i],2) == 0.0:  #special case for nodal plane along x
                xnodes,znodes = numpy.mgrid[-distance_cap:distance_cap:(node_grid*n*1j) , -distance_cap:distance_cap:(node_grid*n*1j)] #floor gives a smaller grid for the nodes
                ynodes = 0*xnodes
            if numpy.round(phi_nodes[i],2) == numpy.round(numpy.pi/2,2):  #special case for nodal plane along y
                ynodes,znodes = numpy.mgrid[-distance_cap:distance_cap:(node_grid*n*1j), -distance_cap:distance_cap:(node_grid*n*1j)] #floor gives a smaller grid for the nodes
                xnodes = 0*ynodes
            xnodes[numpy.where(numpy.abs(xnodes) > distance_cap)] = float('NaN')  #set so the planes stay in the box
            ynodes[numpy.where(numpy.abs(ynodes) > distance_cap)] = float('NaN')  #set so the planes stay in the box
            mlab.mesh(xnodes,ynodes,znodes,representation=Node_mode,color = phi_node_color,opacity=node_opacity,tube_radius=tube_rad)

# Defining a function used for plotting a scalar field for a given orbital #
def Mayavi_Volume(x,y,z,orbital,part):     #orbital is the array containing the orbital information (PhiReal or PhiImag)
    if light_mode == 'black':
        fig = mlab.figure(fgcolor=(1, 1, 1), bgcolor=(0, 0, 0), size=(1200, 800))
    if light_mode == 'white':
        fig = mlab.figure(fgcolor=(0, 0, 0), bgcolor=(1, 1, 1), size=(1200, 800))
    src = mlab.pipeline.scalar_field(x,y,z,orbital)     #generates the data for the 3D density coloring
    volume = mlab.pipeline.volume(src)                  #generates a 3D volume with density coloring defined by src
    if part == 'magnitude':
        Modify_Color_Opacity(volume,orbital,True)  #setting up the custom colormap values and opacities
    if part != 'magnitude':
        Modify_Color_Opacity(volume,orbital,False)  #setting up the custom colormap values and opacities
    if nice_labels:
        axis_ranges = [-(rmax-rmax%10), rmax-rmax%10,-(rmax-rmax%10), rmax-rmax%10,-(rmax-rmax%10), rmax-rmax%10]  #makes axis labels a multiple of 10
    if nice_labels == False:
        axis_ranges = [-rmax, rmax,-rmax, rmax,-rmax, rmax]  #real axes
    axes = mlab.axes(volume, nb_labels=5, xlabel='x (a.u.)', ylabel='y (a.u.)', zlabel='z (a.u.)',line_width = 3, ranges = axis_ranges)
    axes.axes.label_format = '%-#.1f'   #rounds label number to 1 decimal place using f strings
    check = plot_nodes or plot_radial_nodes or plot_theta_nodes or plot_phi_nodes  #checks if the user requested nodal surfaces
    if check and part == 'real':
        Mayavi_Nodes(r_nodes,theta_nodes,real_phi_nodes)  #plots nodal surfaces
    if check and part == 'imag':
        Mayavi_Nodes(r_nodes,theta_nodes,imag_phi_nodes)  #plots nodal surfaces
    if box_outline == True:
        mlab.outline()  #plots a box around the orbital at the edges of the axes
    if part == 'magnitude' and plot_mode == 'original':
        mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=' + str(m) + ' Orbital Magnitude', size=0.2)
    if part == 'real':  #using part to decide which orbital name goes on the plot
        if plot_mode == 'orbital':      
            mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=' + str(m) + ' Orbital', size=0.2)
        if plot_mode == 'probability density':
            mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=' + str(m) + ' Probability Density', size=0.2)
        if plot_mode == 'original':
            mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=' + str(m) + ' Orbital Real Part', size=0.2)
    if part == 'imag':
        if plot_mode == 'orbital':
            mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=-' + str(m) + ' Orbital', size=0.2)
        if plot_mode == 'probability density':
            mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=-' + str(m) + ' Probability Density', size=0.2)
        if plot_mode == 'original':
            mlab.title('n='+ str(n) + ', l=' + str(l) + ', m=' + str(m) + ' Orbital Imaginary Part', size=0.2)
    mlab.view(azimuth=45, elevation = 90, distance = 8*rmax, focalpoint = numpy.array([0.0,0.0,0.0]))  #default viewing angle
print('Python functions successfully defined.')

# Block 6 #
if n-l-1 > 0 and (plot_nodes or plot_radial_nodes):  #Only looks for radial nodes if there are supposed to be some
    r_nodes = scipy.special.roots_genlaguerre(n-l-1,2*l+1, mu=False)[0]*n/2 #easy implemented way to get the roots
    print('The following radial nodes were found (a.u.): ',r_nodes)
    if len(r_nodes) < n-l-1 or len(r_nodes) > n-l-1:
        print('### Warning: ' + str(len(r_nodes)) + ' roots found, ' + str(n-l-1) + ' expected. ###')
if n-l-1 > 0 and (plot_nodes or plot_radial_nodes) == False:
    print('No radial nodes were requested.')
    r_nodes = []
if n-l-1 == 0:
    print('There are no radial nodes since n-l-1 = 0.')
    r_nodes = []
    
if l-abs(m) > 0 and (plot_nodes or plot_theta_nodes):  #Only looks for theta nodes if there are supposed to be some
    theta = sym.symbols('r') #We use sympy to find the locations of the nodes
    Plm_sym = sym.functions.special.polynomials.assoc_legendre(l,m,sym.cos(theta))
    Plm2_sym = Plm_sym**2  #square is easier to search
    solutions = numpy.array([float(sym.nsolve(Plm2_sym,(i+1)*numpy.pi/(3*n),verify=False)) for i in range(3*n)]) #scanning lots of guess values
    solutions = solutions[numpy.where(numpy.round(solutions,1) < numpy.round(numpy.pi,3))] #Throw out false solutions greater than pi
    solutions = solutions[numpy.where(numpy.round(solutions,1) > 0.0)] #Throw out false solutions less than 0
    solutions = numpy.array(list(set(numpy.round(solutions,3)))) #Trick to get rid of redundancy of solutions
    solutions = solutions[numpy.argsort(solutions)]  #sorting the solutions
    solutions = solutions[numpy.where(numpy.round(solutions,2) != 0.0)]  #used to check if the solutions are actually solutions
    solutions = solutions[numpy.where(numpy.round(solutions,2) != numpy.round(numpy.pi,2))]  #used to check if the solutions are actually solutions
    success = numpy.round([float(Plm2_sym.subs(theta,solutions[i])) for i in range(len(solutions))],1) #checking which solutions work
    theta_nodes = solutions[numpy.where(numpy.round(success,2) == 0.0)]   #unique solutions that work are the theta nodes
    print('the following theta nodes were found (units of pi): ',theta_nodes/numpy.pi)
    if len(theta_nodes) < l-abs(m) or len(theta_nodes) > l-abs(m):
        print('### Warning: ' + str(len(theta_nodes)) + ' roots found, ' + str(l-m) + ' expected. ###')
if l-abs(m) > 0 and (plot_nodes or plot_theta_nodes) == False:
    print('No theta nodes were requested.')
    theta_nodes = []
if l-abs(m) == 0:
    print('There are no theta nodes since l-m = 0.')
    theta_nodes = []

if abs(m) > 0 and (plot_nodes or plot_phi_nodes):  #Only calculated phi nodes if there are supposed to be some
    phi_step = numpy.pi/abs(m)  #steps between nodes along phi
    real_phi_nodes = numpy.array([i*phi_step + phi_step/2 for i in range(abs(m))])
    imag_phi_nodes = numpy.array([i*phi_step for i in range(abs(m))]) 
    print('The real part phi nodes are (units of pi): ', real_phi_nodes/numpy.pi)
    print('The imaginary part phi nodes are (units of pi): ', imag_phi_nodes/numpy.pi)
if abs(m) > 0 and (plot_nodes or plot_phi_nodes) == False:
    print('No phi nodes were requested.')
    real_phi_nodes = []
    imag_phi_nodes = []
if m == 0:
    print('There are no phi nodes since m = 0.')
    real_phi_nodes = []
    imag_phi_nodes = []

# Block 7 #
if wavefunction_mode != 'original' and m < 0:
    print('m has been changed to a positive value to match phi_n,l,m definition.')
    m = abs(m)      #keeping m positive since we don't need negative values for phi

# Generating the Colormap Data #
if reverse_colormap == True:
    COL = MplColorHelper(colormap + '_r', 0, 1)   #_r reverses the order of the colormap, which can look better
if reverse_colormap == False:
    COL = MplColorHelper(colormap, 0, 1)
color_array = numpy.zeros((colormap_points,3)) #RGB values used for the hacky colormap; must be odd rows to catch zero! ~11 gives smooth color gradients
COL = COL.get_rgb(numpy.linspace(0,1,len(color_array)))  #turning COL into the actual colormap RGB values
color_array = COL[:,0:3]    #filling array of RGB data for the custom colormap
reverse_color_array = color_array[::-1]
if reverse_colormap:
    print('Colormap data for ' + str(colormap) + '_r generated.')
else:
    print('Colormap data for ' + str(colormap) + ' generated.')

# Generating the Coordinate Grid Data #
x,y,z = numpy.mgrid[-rmax:rmax:(N*1j) , -rmax:rmax:(N*1j) , -rmax:rmax:(N*1j)]  #setting up N*N*N mgrid arrays for x,y,z with bounds of +/-rmax
theta = x*0  #initializing array to be the right shape
phi = x*0
r = numpy.sqrt(x**2 + y**2 + z**2)                      #making the radial values
R = numpy.sqrt(x**2 + y**2)                             #useful for making phi values (distance in x-y plane)
indices = numpy.where(r>0)                              #selecting indices which will not break theta
theta[indices] = numpy.arccos(z[indices]/r[indices])    #making the theta values that are allowed; 0s elsewhere
indices = numpy.where(R>0)                              #selecting indices which will not break phi
phi[indices] = numpy.arccos(x[indices]/R[indices])      #making the phi grid; still need to fix values for y < 0
phi[numpy.where(y<0)] = -phi[numpy.where(y<0)]          #multiplying phi value by -1 if y < 0 (necessary; will not work if you don't use mgrid!)
print('Coordinate grids for x,y,z and r,theta,phi generated.')

# Generating the Function Data #
Cnl = math.sqrt((2/n)**3*math.factorial(n-l-1)/(2*n*math.factorial(n+l)))       #normalization constant for Rnl
if wavefunction_mode != 'original':
    Nlm = math.sqrt((2*l+1)*math.factorial(l-abs(m))/(4*numpy.pi*math.factorial(l+abs(m))))   #normalization constant for Yml
if wavefunction_mode == 'original':
    Nlm = math.sqrt((2*l+1)*math.factorial(l-m)/(4*numpy.pi*math.factorial(l+m)))   #normalization constant for Yml
rho = 2*r/n         #shortcut saving calculations
LaP = scipy.special.genlaguerre(n-l-1,2*l+1)(rho)       #radial wavefunction's Laguerre polynomial
Rnl = Cnl*numpy.exp(-rho/2)*rho**l*LaP                  #radial component of wavefunction
LeP = (-1)**m*scipy.special.lpmv(abs(m),l,numpy.cos(theta))  #Legendre polynomial for angular component of wavefunction; Condon-Shortley phase is removed
Yml = Nlm*LeP*numpy.exp(complex(0,1)*m*phi)             #angular component of wavefunction is a spherical harmonic function
print('Component functions generated on the grid.')

# Block 8 #
# Generating the real part orbital phi_n,l,m or the real part of psi_n,l,m on the grid #
if plot_mode != 'original':
    print('Generating Real Part Orbital phi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
if plot_mode == 'original':
    print('Generating Real Part of Orbital psi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
PhiReal = x*0   #making empty array that has the right shape for real part orbital
PhiMax = numpy.amax(numpy.abs(Rnl*Yml))     #use the max value to scale the wavefunction for simplicity
indices = numpy.where(numpy.abs(Rnl*numpy.real(Yml))/PhiMax <= cutoff) #indices where the wavefunction is below the cutoff
PhiReal[indices] = 0.0             #turning off values for real part orbital phi_n,l,m
indices = numpy.where(numpy.abs(Rnl*numpy.real(Yml))/PhiMax > cutoff)  #indices that we will actually use to generate the orbital
if wavefunction_mode == 'scaled':
    PhiReal[indices] = Rnl[indices]*numpy.real(Yml[indices])/PhiMax  #real part orbital phi_n,l,m scaled by the max value
if wavefunction_mode == 'actual':
    if plot_mode == 'original':
        PhiReal[indices] = Rnl[indices]*numpy.real(Yml[indices])     #real part orbital psi_n,l,m; called Phi to save code
    if plot_mode == 'orbital':
        PhiReal[indices] = math.sqrt(2)*Rnl[indices]*numpy.real(Yml[indices])     #real part orbital phi_n,l,m
if (m == 0 and m0_cut):  #turning off the first quadrant for s orbitals or orbitals with m = 0
    indices = numpy.where(numpy.logical_and(x>0,y>0))
    PhiReal[indices] = 0.0
if plot_mode != 'original':
    print('...Real Part Orbital phi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated on grid.')
    print('Generating 3-D Volume Data for phi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
if plot_mode == 'original':
    print('...Real Part of Orbital psi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated on grid.')
    print('Generating 3-D Volume Data for psi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
if plot_mode == 'orbital' or plot_mode == 'original':
    Mayavi_Volume(x,y,z,PhiReal,'real')    #creating the 3-D volumetric plot for the real part orbital
if plot_mode == 'probability density':
    Prob = x*0  #Initializing array
    ProbMax = numpy.amax(Rnl**2*numpy.real(Yml)**2)    #highest value of the probability density
    indices = numpy.where(numpy.abs(Rnl**2*numpy.real(Yml)**2)/ProbMax <= cutoff) 
    Prob[indices] = 0.0    #cutting off probability density function
    indices = numpy.where(Rnl**2*numpy.real(Yml)**2/ProbMax > cutoff)
    #sign = numpy.sign(numpy.real(Yml[indices]))*numpy.sign(Rnl[indices])  #phases from the orbital
    if wavefunction_mode == 'scaled':
        #Prob[indices] = sign*Rnl[indices]**2*numpy.real(Yml[indices])**2/ProbMax  #retains the phase information
        Prob[indices] = Rnl[indices]**2*numpy.real(Yml[indices])**2/ProbMax
    if wavefunction_mode == 'actual':
        #Prob[indices] = sign*Rnl[indices]**2*numpy.real(Yml[indices])**2  #retains the phase information
        Prob[indices] = Rnl[indices]**2*numpy.real(Yml[indices])**2
    if (m == 0 and m0_cut):  #turning off the first quadrant for s orbitals or orbitals with m = 0
        indices = numpy.where(numpy.logical_and(x>0,y>0))
        Prob[indices] = 0.0
    Mayavi_Volume(x,y,z,Prob,'real')    #creating the 3-D volumetric plot for the real part orbital probabilty density
if plot_mode != 'original':
    print('...3-D Volume Data for phi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated.')
if plot_mode == 'original':
    print('...3-D Volume Data for psi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated.')
    
# Generating the imaginary part orbital phi_n,l,-m or the imaginary part of psi_n,l,m on the grid #
if abs(m) > 0:   #no need to plot imag orbital if m = 0
    PhiImag = x*0   #making 0 array that has the right shape for imag part orbital
    # Generating the imaginary part orbital on the grid #
    if plot_mode != 'original':
        print('Generating Imaginary Part Orbital phi_' + str(n) + ',' + str(l) + ',-' + str(m) + '...')
    if plot_mode == 'original':
        print('Generating Imaginary Part of Orbital psi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
    indices = numpy.where(numpy.abs(Rnl*numpy.imag(Yml))/PhiMax >= cutoff)  #only populates values greater than cutoff
    if wavefunction_mode == 'scaled':
        PhiImag[indices] = Rnl[indices]*numpy.imag(Yml[indices])/PhiMax   #imag part orbital phi_n,l,-m scaled by the max value
    if wavefunction_mode == 'actual':
        if plot_mode == 'original':
            PhiImag[indices] = Rnl[indices]*numpy.imag(Yml[indices])   #imag part orbital phi_n,l,-m
        if plot_mode == 'orbital':
            PhiImag[indices] = math.sqrt(2)*Rnl[indices]*numpy.imag(Yml[indices])   #imag part orbital phi_n,l,-m
    if plot_mode != 'original':
        print('...Imaginary Part Orbital phi_' + str(n) + ',' + str(l) + ',-' + str(m) + ' generated on grid.')
        print('Generating 3-D Volume Data for phi_' + str(n) + ',' + str(l) + ',-' + str(m) + '...')
    if plot_mode == 'original':
        print('...Imaginary Part of Orbital psi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated on grid.')
        print('Generating 3-D Volume Data for psi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
    if plot_mode == 'orbital' or plot_mode == 'original':
        Mayavi_Volume(x,y,z,PhiImag,'imag')     #creating the 3-D volumetric plot for the imaginary part orbital
    if plot_mode == 'probability density':
        Prob = x*0  #Initializing array
        ProbMax = numpy.amax(Rnl**2*numpy.imag(Yml)**2)    #highest value of the probability density
        indices = numpy.where(numpy.abs(Rnl**2*numpy.imag(Yml)**2)/ProbMax <= cutoff) 
        Prob[indices] = 0.0    #cutting off probability density function
        indices = numpy.where(Rnl**2*numpy.imag(Yml)**2/ProbMax > cutoff)
        #sign = numpy.sign(numpy.imag(Yml[indices]))*numpy.sign(Rnl[indices])  #phases from the orbital
        if wavefunction_mode == 'scaled':
            #Prob[indices] = sign*Rnl[indices]**2*numpy.imag(Yml[indices])**2/ProbMax  #retains the phase information
            Prob[indices] = sign*Rnl[indices]**2*numpy.imag(Yml[indices])**2/ProbMax
        if wavefunction_mode == 'actual':
            #Prob[indices] = sign*Rnl[indices]**2*numpy.imag(Yml[indices])**2  #retains the phase information
            Prob[indices] = sign*Rnl[indices]**2*numpy.imag(Yml[indices])**2
        if (m == 0 and m0_cut):  #turning off the first quadrant for s orbitals or orbitals with m = 0
            indices = numpy.where(numpy.logical_and(x>0,y>0))
            Prob[indices] = 0.0
        Mayavi_Volume(x,y,z,Prob,'imag')    #creating the 3-D volumetric plot for the real part orbital probabilty density
    if plot_mode != 'original':
        print('...3-D Volume Data for phi_' + str(n) + ',' + str(l) + ',-' + str(m) + ' generated.')
    if plot_mode == 'original':
        print('...3-D Volume Data for psi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated.')

if plot_mode == 'original' and abs(m) > 0:
    print('Generating 3-D Volume Data for Magnitude of psi_' + str(n) + ',' + str(l) + ',' + str(m) + '...')
    Mayavi_Volume(x,y,z,numpy.abs(PhiReal + 1j*PhiImag),'magnitude')  #extra plot for the magnitude of the orbital
    print('...3-D Volume Data for Magnitude of psi_' + str(n) + ',' + str(l) + ',' + str(m) + ' generated.')
        
mlab.show()     #only once the program gets here does it actually show anything; this can only be called once or things break!

# Block 9 #
# Generating plot of the radial wavefunction and radial probability density #
plt.rcParams['figure.dpi'] = 200    #makes the figure higher resolution than the Jupyter default
try:
    plt.rcParams['text.usetex'] = True  #allows use of LaTeX in labels
except:
    print('No LaTeX Distribution Detected, Using Plain Text Labels.')
plt.figure(num=1, figsize=(5, 3))   #size is specified because Jupyter randomly resizes things otherwise
r = numpy.linspace(0,rmax,2*N)
rho = 2*r/n   #convenient for calculations
LaP = scipy.special.genlaguerre(n-l-1,2*l+1)(rho)   #radial wavefunction Laguerre polynomial part
Rnl = Cnl*numpy.exp(-rho/2)*rho**l*LaP              #radial component of wavefunction
Pnl = r**2*Rnl**2  #radial probability density
if wavefunction_mode == 'scaled':
    Rnl = Rnl/max(Rnl)   #Rescaling to max value 1
    Pnl = Pnl/max(Pnl)   #Rescaling to max value 1
plt.plot([0,rmax],[0,0], color = 'k', linestyle = '--')
try:
    plt.plot(r, Rnl, lw=2, color = 'blue', label = r'$R_{' + str(n) + ',' + str(l) + '}(r)$') 
except:
    plt.plot(r, Rnl, lw=2, color = 'blue', label = 'Radial Part') 
if (plot_nodes or plot_radial_nodes) or plot_mode == 'probability':
    try:
        plt.plot(r, Pnl,  lw=2, color = 'red', label = r'$P_{' + str(n) + ',' + str(l) + '}(r)$')
    except:
        plt.plot(r, Pnl,  lw=2, color = 'red', label = 'radial part')
if (plot_nodes or plot_radial_nodes) and len(r_nodes) > 0:
    plt.scatter(r_nodes,numpy.zeros(len(r_nodes)),marker='o',color = 'green',facecolors='none')
plt.xlim([0,rmax])
try:
    plt.xlabel(r'$r$ (a.u.)')
except:
    plt.xlabel('r (a.u.)')
if wavefunction_mode == 'scaled':
    plt.ylim([min(Rnl)-0.1,1.1])
    if (plot_nodes or plot_radial_nodes):
        try:
            plt.ylabel(r'$R_{' + str(n) + ',' + str(l) + '}(r), P_{' + str(n) + ',' + str(l) + '}(r)$ (Scaled)')
        except:
            plt.ylabel('Radial Part, Radial Probability (Scaled)')
    else:
        try:
            plt.ylabel(r'$R_{' + str(n) + ',' + str(l) + '}(r)$ (Scaled)')
        except:
            plt.ylabel('Radial Part (Scaled)')
if wavefunction_mode == 'actual':
    plt.ylim([min(Rnl)-0.1*max(max(Rnl),max(Pnl)),1.1*max(max(Rnl),max(Pnl))])
    if (plot_nodes or plot_radial_nodes):
        try:
            plt.ylabel(r'$R_{' + str(n) + ',' + str(l) + '}(r), P_{' + str(n) + ',' + str(l) + '}(r)$')
        except:
            plt.ylabel('Radial Part, Radial Probability')
    else:
        try:
            plt.ylabel(r'$R_{' + str(n) + ',' + str(l) + '}(r)$')
        except:
            plt.ylabel('Radial Part')
plt.legend(bbox_to_anchor=(1.05, 0.6), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.show()  #must not be shown until after the mlab.show() or else the program will freeze

# Block 10 #
# Generating plot of the theta angular wavefunction #
if l == 0:
    print('Theta part is trivial since l = 0 gives P_0^0 = 1')
if l > 0:
    plt.figure(num=2, figsize=(5, 3))
    theta = numpy.linspace(0,numpy.pi,N)
    LeP = (-1)**m*Nlm*scipy.special.lpmv(m,l,numpy.cos(theta))  #We give the normalization constant to theta part
    if wavefunction_mode == 'scaled':
        LeP = LeP/max(LeP)
    plt.plot([0,1],[0,0], color = 'k', linestyle = '--')
    try:
        plt.plot(theta/numpy.pi, LeP, lw=2, color = 'blue', label = r'$P_' + str(l) + '^' + str(m) + r'(\cos(\theta))$') 
    except:
        plt.plot(theta/numpy.pi, LeP, lw=2, color = 'blue', label = 'theta Part') 
    if (plot_nodes or plot_theta_nodes) or plot_mode == 'probability': 
        try:
            plt.plot(theta/numpy.pi, LeP**2, lw=2, color = 'red', label = r'$P_' + str(l) + '^' + str(m) + r'(\cos(\theta))^2$')
        except:
            plt.plot(theta/numpy.pi, LeP**2, lw=2, color = 'red', label = 'Square of theta part')
    if (plot_nodes or plot_theta_nodes) and len(theta_nodes) > 0:
        plt.scatter(theta_nodes/numpy.pi,numpy.zeros(len(theta_nodes)),marker='o',color = 'green',facecolors='none')
    plt.xlim([0,1])
    try:
        plt.xlabel(r'$\theta$ ($\pi$)')
    except:
        plt.xlabel('theta (pi)')
    if wavefunction_mode == 'scaled':
        plt.ylim([min(LeP)-0.1*max(LeP),1.1*max(LeP)])
        if (plot_nodes or plot_theta_nodes):
            try:
                plt.ylabel(r'$P_' + str(l) + '^' + str(m) + r'(\cos(\theta))$, $P_' + str(l) + '^' + str(m) + r'(\cos(\theta))^2$ (Scaled)')
            except:
                plt.ylabel('theta Part, Square of theta Part (Scaled)')
        else:
            try:
                plt.ylabel(r'$P_' + str(l) + '^' + str(m) + r'(\cos(\theta))$ (Scaled)')
            except:
                plt.ylabel('theta Part (Scaled)')
    if wavefunction_mode == 'actual':
        plt.ylim([min(LeP)-0.1*max(max(LeP),max(LeP**2)),1.1*max(max(LeP),max(LeP**2))])
        if (plot_nodes or plot_theta_nodes):
            try:
                plt.ylabel(r'$P_' + str(l) + '^' + str(m) + r'(\cos(\theta))$, $P_' + str(l) + '^' + str(m) + r'(\cos(\theta))^2$')
            except:
                plt.ylabel('theta Part, Square of theta Part')
        else:
            try:
                plt.ylabel(r'$P_' + str(l) + '^' + str(m) + r'(\cos(\theta))$')
            except:
                plt.ylabel('theta Part')
    plt.legend(bbox_to_anchor=(1.05, 0.6), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.show()  #must not be shown until after the mlab.show() or else the program will freeze

# Block 11 #
# Generating plot of the phi angular wavefunction #
if m == 0:
    print('Phi part is trivial since m = 0 gives exp(0) = 1')
if abs(m) > 0 :
    plt.figure(num=3, figsize=(5,3))
    phi = numpy.linspace(0,2*numpy.pi,2*N)
    F = numpy.exp(complex(0,1)*m*phi)  #the wavefunction for phi
    plt.plot([0,2],[0,0], color = 'k', linestyle = '--')
    if plot_nodes == True:
        plt.plot([1,1],[-1.3,1.3], color = 'green', linestyle = '--')
    if (plot_nodes or plot_phi_nodes):
        if len(real_phi_nodes) > 0:
            plt.scatter(real_phi_nodes/numpy.pi,numpy.zeros(len(real_phi_nodes)),marker='o',color = 'green',facecolors='none')
        if len(imag_phi_nodes) > 0:
            plt.scatter(imag_phi_nodes/numpy.pi,numpy.zeros(len(imag_phi_nodes)),marker='s',color = 'green',facecolors='none')
    if m > 1:
        try:
            plt.plot(phi/numpy.pi, numpy.real(F), lw=2, color = 'blue', label = r'Re$(e^{' + str(m) + 'i \phi})$')
            plt.plot(phi/numpy.pi, numpy.imag(F), lw=2, color = 'red', label = r'Im$(e^{' + str(m) + 'i \phi})$')
        except:
            plt.plot(phi/numpy.pi, numpy.real(F), lw=2, color = 'blue', label = 'Real Part for phi')
            plt.plot(phi/numpy.pi, numpy.imag(F), lw=2, color = 'red', label = 'Imaginary Part for phi')
    if m == 1:
        try:
            plt.plot(phi/numpy.pi, numpy.real(F), lw=2, color = 'blue', label = r'Re$(e^{i \phi})$')
            plt.plot(phi/numpy.pi, numpy.imag(F), lw=2, color = 'red', label = r'Im$(e^{i \phi})$')
        except:
            plt.plot(phi/numpy.pi, numpy.real(F), lw=2, color = 'blue', label = 'Real Part for phi')
            plt.plot(phi/numpy.pi, numpy.imag(F), lw=2, color = 'red', label = 'Imaginary Part for phi')
    plt.xlim([0,2])
    plt.ylim([-1.3,1.3])
    try:
        plt.xlabel(r'$\phi$ ($\pi$)')
    except:
        plt.xlabel('phi (pi)')
    if m > 1:
        try:
            plt.ylabel(r'Re$(e^{' + str(m) + 'i \phi})$,Im$(e^{i' + str(m) + '\phi})$')
        except:
            plt.ylabel('Real and Imaginary Parts for phi')
    if  m == 1:
        try:
            plt.ylabel(r'Re$(e^{i \phi})$,Im$(e^{i \phi})$')
        except:
            plt.ylabel('Real and Imaginary Parts for phi')
    plt.legend(bbox_to_anchor=(1.05, 0.6), loc=2, borderaxespad=0.)
    plt.tight_layout()
    plt.show()  #must not be shown until after the mlab.show() or else the program will freeze

# Block 12 #
plt.figure(num=4, figsize=(5, 3))
if wavefunction_mode == 'scaled':
    values = numpy.linspace(-1,1,colormap_points) 
if wavefunction_mode == 'actual':
    if plot_mode == 'orbital' or plot_mode == 'original':
        values = numpy.linspace(-PhiMax,PhiMax,colormap_points)
    if plot_mode == 'probability density':
        values = numpy.linspace(-PhiMax**2,PhiMax**2,colormap_points)
plt.plot([min(values),max(values)],[0,0], color = 'k', linestyle = '--')
opacities = [opacity_factor*((numpy.abs(-colormap_points//2+i+1)/(colormap_points//2))**opacity_exponent)+opacity_shift for i in range(colormap_points)]
plt.scatter(values,color_array[:,0], lw=1, color = 'red', label = 'Red values for ' + colormap,marker='d',facecolors='none')
plt.scatter(values,color_array[:,1], lw=1, color = 'green', label = 'Green values for ' + colormap,marker='s',facecolors='none')
plt.scatter(values,color_array[:,2], lw=1, color = 'blue', label = 'Blue values for ' + colormap,marker='*',facecolors='none')
plt.scatter(values,opacities, lw=1, color = 'black', label = 'Opacity values',marker='o',facecolors='none')
plt.plot(values,color_array[:,0], lw=1, color = 'red')
plt.plot(values,color_array[:,1], lw=1, color = 'green')
plt.plot(values,color_array[:,2], lw=1, color = 'blue')
plt.plot(values,opacities, lw=1, color = 'black')
if plot_mode == 'probability density':
    plt.plot([1,1],[-1.1,1.1], color = 'black', linestyle = '--')  #shows that probability only uses positive values
plt.xlim([min(values),max(values)])
plt.ylim([-0.1,1.1])
if plot_mode == 'orbital':
    label_string = 'Orbital Value'
if plot_mode == 'probability density':
    label_string = 'Probability Density Value'
if plot_mode == 'original':
    label_string = 'Orbital Value'
if wavefunction_mode == 'scaled':
    label_string += ' (Scaled)'
plt.xlabel(label_string)
plt.ylabel('Map Value')
plt.legend(bbox_to_anchor=(1.05, 0.7), loc=2, borderaxespad=0.)
plt.tight_layout()
plt.show()  #must not be shown until after the mlab.show() or else the program will freeze
