import os
import numpy as np
import pylorentz

class unit:
    m=100
    cm=1
    

def frame_transform(four_vector_in, boost_direction, beta=None, gamma=None):
    """
    Transform the four-vector (e, p_x, p_y, p_z,) with boost_direction and beta/gamma.
    
    INPUT:
    ---
    e, p_x, p_y, p_z: 
        four-vector to be transformed. Can be momentum or what ever 4-vector
    boost_direction:
        three-vector, the relative velocity of the frame to be transformed to
    beta/gamma:
        at least one needs to be given
        
    RETURN
    ---
    four_vector_out:
        the transformed four-vector
    """
    x, y, z = boost_direction
    e, p_x, p_y, p_z = four_vector_in
    four_vector_out = pylorentz.Momentum4(e, p_x, p_y, p_z).boost(x, y, z, beta=beta, gamma=gamma)
    return four_vector_out
    
    
def twobody_decay(mass, three_momentum, decayproduct_pid=[13,13], rand_seed=None):
    """
    INPUT
    ---
    mass: 
        float, mass of the primary particle [MeV]
    three_momentum: 
        list, three momentum of the primary partile [MeV]
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    """
    mass_map={13: 105.7,    # mu-
            11:0.511,       # e-
              211: 139.57,  # pi+/-
              -211: 139.57, # pi+/-
             }
    
    m1=mass_map[decayproduct_pid[0]]
    m2=mass_map[decayproduct_pid[1]]
    M=mass
    P=np.sqrt(np.sum(np.power(three_momentum,2)))
    
    # First, 0->1+2 decay in the rest frame of 0:
    E1 = (M**2+m1**2-m2**2)/(2*M)
    E2 = (M**2+m2**2-m1**2)/(2*M)
    p1 = ((M**2+m1**2-m2**2)**2- 4*M**2*m1**2)**0.5/(2*M)
    p2 = ((M**2+m2**2-m1**2)**2- 4*M**2*m2**2)**0.5/(2*M)
    
    # Generate a randomized direction
    rng = np.random.default_rng(seed=rand_seed)
    p1_x = rng.normal(0,1)
    p1_y = rng.normal(0,1)
    p1_z = rng.normal(0,1)
    mag = np.sqrt(p1_x**2+p1_y**2+p1_z**2)
    p1_x = p1_x/mag*p1
    p1_y = p1_y/mag*p1
    p1_z = p1_z/mag*p1
    
    # Boost the 4-momentum
    boost_direction = -np.array(three_momentum)
    boost_gamma = np.sqrt(M**2+P**2)/M
    vec1 = frame_transform([E1, p1_x, p1_y, p1_z], boost_direction, gamma=boost_gamma)
    vec2 = frame_transform([E1, -p1_x, -p1_y, -p1_z], boost_direction, gamma=boost_gamma)
    
    return vec1, vec2

def twobody_decay_MA(mass, abs_momentum, vertex_xyz, decayproduct_pid=[13,13], rand_seed=None):
    """
    INPUT
    ---
    mass: 
        float, mass of the primary particle [MeV]
    abs_momentum: 
        float, absolute momentum of the primary partile [MeV]
    vertex_xyz: in uni of [cm]
        [x,y,z] in **CMS** coordinates! It's the same coordinates as the HIT_x,y,z in simulation output.
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    
    RETURN
    ---
    p4vec_1,p4vec_2:
        Momentum four-vector of the two decay products
    
    """
    r = np.sqrt(np.sum(np.power(vertex_xyz,2)))
    px = abs_momentum/r*vertex_xyz[0]
    py = abs_momentum/r*vertex_xyz[1]
    pz = abs_momentum/r*vertex_xyz[2]
    
    p4vec_1,p4vec_2 = twobody_decay(mass, [px,py,pz], decayproduct_pid=decayproduct_pid, rand_seed=rand_seed) 
    
    return p4vec_1,p4vec_2

def multibody_decay_MA(mass, abs_momentum, vertex_xyz, decayproduct_p4vec_list):
    boost_direction = -np.array(three_momentum)
    boost_gamma = np.sqrt(mass**2+abs_momentum**2)/mass
    
    p4vec_transformed=np.array([frame_transform(particle, boost_direction, gamma=boost_gamma) for particle in decayproduct_p4vec_list])
    return p4vec_transformed 

def generate_sim_script_filereader(events_properties_filename, script_path=None):
    """
    Create simulation script for filereader generator based on the events database file.
    """
    if script_path is None:
        script_path = os.path.splitext(events_properties_filename)[0]+".mac"
        
    with open(events_properties_filename, 'r') as file:
        nlines_read=0
        found_nevents=False
        while (not found_nevents) and nlines_read<=10:
            line_content = file.readline().split()
            nlines_read+=1
            if "nevents" in line_content:
                try:
                    nevents = int(line_content[-1])
                    found_nevents = True
                    break
                except:
                    print("Could not read number of events")
            

    script = "/det/select Box \n"
    script+= "/gen/select file_reader \n"
    script+= f"/gen/file_reader/pathname {events_properties_filename}\n"
    script+= f"/run/beamOn {nevents}"

    with open(script_path, 'w') as file:
        file.write(script)
        
    print("Script saved at", script_path)
    
    return script_path

def generate_twobody_decay_file(filename, mass, abs_momentum, vertex_xyz, decayproduct_pid=[13,13], rand_seed=None, Nevents=10000, OVER_WRITE=False, which_coordinate = "CMS"):
    """
    INPUT
    ---
    filename:
        str, the filename of the event description file to be generated
    mass: 
        float, mass of the primary particle [MeV]
    abs_momentum: 
        float, absolute momentum of the primary partile [MeV]
    vertex_xyz: [x,y,z] or [[x0,x1], [y0,y1], [z0,z1]], in unit of [cm]
        [x,y,z] in **CMS** coordinates! It's the same coordinates as the HIT_x,y,z in simulation output.
        [[x0,x1], [y0,y1], [z0,z1]]: ranges of xyz, will generate uniform distributed xyz. 
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    which_coordinate:
        "CMS" or "detector"
    """    
    
    if os.path.exists(filename):
        print("File exists!")
        if not OVER_WRITE:
            print("Please change filename, or set OVER_WRITE=True")
            return
        
    # Transform vertex position to detector coordinate to feed to simulation
    vertex_xyz = np.array(vertex_xyz)
    if vertex_xyz.ndim==1:
        if which_coordinate=="CMS":
            vertex_xyz_det = np.array([vertex_xyz[2],-vertex_xyz[0],-vertex_xyz[1]+85.47*unit.m])*10 # turn into mm
            vertex_xyz_cms = vertex_xyz*10
        else:
            vertex_xyz_det = vertex_xyz*10
            vertex_xyz_cms = np.array([-vertex_xyz[1],-vertex_xyz[2]+85.47*unit.m,vertex_xyz[0]])*10
        
    
    rng = np.random.default_rng(seed=rand_seed)
    
    
    with open(filename, "w") as file:
        # first, write the total number of events
        primary_particle_PID = -1
        r = np.sqrt(np.sum(np.power(vertex_xyz_cms,2)))
        primary_particle_px = abs_momentum/r*vertex_xyz_cms[0]
        primary_particle_py = abs_momentum/r*vertex_xyz_cms[1]
        primary_particle_pz = abs_momentum/r*vertex_xyz_cms[2]     
        
        file.write(f"# nevents {Nevents}\n\n")
        for i in range(Nevents):
            file.write(f"n {i}  \t {primary_particle_PID}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {primary_particle_px}\t {primary_particle_py}\t {primary_particle_pz}\n")
            rand_seed_i = None if rand_seed is None else rand_seed+i
            
            # Generate a random vertex if a range was given
            if vertex_xyz.ndim==2:
                vertex_temp = np.array([rng.uniform(*vertex_xyz[0]), rng.uniform(*vertex_xyz[1]), rng.uniform(*vertex_xyz[2])])
                if which_coordinate=="CMS":
                    # Convert to dtector coordinate and # turn into mm
                    vertex_xyz_det = np.array([vertex_temp[2],-vertex_temp[0],-vertex_temp[1]+85.47*unit.m])*10 
                    vertex_xyz_cms = vertex_temp*10
                    
                else:
                    vertex_xyz_det = vertex_temp*10
                    vertex_xyz_cms = np.array([-vertex_temp[1],-vertex_temp[2]+85.47*unit.m,vertex_temp[0]])*10
                    

            p4vec_1,p4vec_2 = twobody_decay_MA(mass, abs_momentum, vertex_xyz_cms, decayproduct_pid=decayproduct_pid, rand_seed=rand_seed_i)
            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {p4vec_1[0]}\t {p4vec_1[1]}\t {-p4vec_1[2]}\n")
            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {p4vec_2[0]}\t {p4vec_2[1]}\t {-p4vec_2[2]}\n")
            
            
            
            
def generate_twobody_vertex_range(filename, abs_momentum, vertex_xyz, theta_range, p_unit_pre=[1,1,-1], phi_range=[0,2*np.pi], decayproduct_pid=[13,13], rand_seed=None, Nevents=10000, OVER_WRITE=False):
    """
    INPUT
    ---
    filename:
        str, the filename of the event description file to be generated
    abs_momentum: 
        float, absolute momentum of the two particles [MeV]
    vertex_xyz: [x,y,z] or [[x0,x1], [y0,y1], [z0,z1]], in unit of [cm]
        [x,y,z] in **DETECTOR** coordinates! It's the same coordinates as the x,y,z of the 'Range' generator in simulation output.
        [[x0,x1], [y0,y1], [z0,z1]]: ranges of xyz, will generate uniform distributed xyz. 
    theta_range:
        [theta_low, theta_high]
    decayproduct_pid:
        [pid1, pid2], PID of the two decay products. Used to decide the mass of the decay particle.
    rand_seed: 
        If None, then fresh, unpredictable entropy will be pulled from the OS. 
    which_coordinate:
        "CMS" or "detector"
    """    
    
    if os.path.exists(filename):
        print("File exists!")
        if not OVER_WRITE:
            print("Please change filename, or set OVER_WRITE=True")
            return
        
    # Transform vertex position to detector coordinate to feed to simulation
    vertex_xyz = np.array(vertex_xyz)
    vertex_xyz_det = vertex_xyz*10
    vertex_xyz_cms = np.array([-vertex_xyz[1],-vertex_xyz[2]+85.47*unit.m,vertex_xyz[0]])*10
        
    rng = np.random.default_rng(seed=rand_seed)
    with open(filename, "w") as file:
        # first, write the total number of events
        file.write(f"# nevents {Nevents}\n\n")
        
        primary_particle_PID = -1
        r = np.sqrt(np.sum(np.power(vertex_xyz_cms,2)))
        primary_particle_px = abs_momentum/r*vertex_xyz_cms[0]
        primary_particle_py = abs_momentum/r*vertex_xyz_cms[1]
        primary_particle_pz = abs_momentum/r*vertex_xyz_cms[2]     
        
        
        for i in range(Nevents):
            file.write(f"n {i}  \t {primary_particle_PID}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {primary_particle_px}\t {primary_particle_py}\t {primary_particle_pz}\n")
            
            rand_seed_i = None if rand_seed is None else rand_seed+i
            
            pvecs=[]
            for ivec in range(2):
                theta = rng.uniform(*theta_range)
                phi = rng.uniform(*phi_range)
                p_unit=np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), -np.cos(theta)])
                p_unit=p_unit/np.sqrt(np.sum(p_unit**2))*p_unit_pre
                p = abs_momentum*p_unit         
                pvecs.append(p)
            
            # Generate a random vertex if a range was given
            if vertex_xyz.ndim==2:
                vertex_temp = np.array([rng.uniform(*vertex_xyz[0]), rng.uniform(*vertex_xyz[1]), rng.uniform(*vertex_xyz[2])])
                vertex_xyz_det = vertex_temp*10
                    

            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {pvecs[0][0]}\t {pvecs[0][1]}\t {-pvecs[0][2]}\n")
            file.write(f"\t {decayproduct_pid[0]}\t  {vertex_xyz_det[0]}\t  {vertex_xyz_det[1]}\t {vertex_xyz_det[2]}\t {pvecs[1][0]}\t {pvecs[1][1]}\t {-pvecs[1][2]}\n")            
            
            
            