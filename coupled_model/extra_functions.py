import numpy as np
import matplotlib.pyplot as plt

def plot(coupling, meshnr, dt, T, k6, bigE, str_C1, bcs):

    #temp_der_stats = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.npy')
    #temp_stats_minmax = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-minmax-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.npy')
    #temp_stats = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.npy')
    #variables = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-oth-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.npy')
    temp_der_stats = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')
    temp_stats_minmax = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-minmax-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')
    temp_stats = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')
    variables = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'-oth-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')


    # plot the temporal statistics
    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten() 
    axs[0].plot(temp_der_stats[4,:], temp_der_stats[0,:])
    axs[1].plot(temp_der_stats[4,:], temp_der_stats[1,:])
    axs[2].plot(temp_der_stats[4,:], temp_der_stats[2,:])
    #axs[3].plot(temp_der_stats[4,:], temp_der_stats[3,:])
    axs[3].plot(temp_der_stats[4,1:], temp_der_stats[3,1:])
    for i in range(4):
        #axs[i].set_yscale('log')
        axs[i].set_xlabel('Time')
        axs[i].set_ylabel('L2-norm')
        axs[i].grid(alpha=0.3, which='both', linestyle=':')
    axs[0].set_title('$||\partial_t c_a||_{L^2(Y)}$')
    axs[1].set_title('$||\partial_t c_d||_{L^2(Y)}$')
    axs[2].set_title('$||\partial_t p_a||_{L^2(\Gamma)}$')
    axs[3].set_title('$||\partial_t u||_{L^2(Y)}$')
    #plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.png')
    plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.png')
    #plt.show()
    plt.close()

    # plot the temporal statistics
    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten()
    axs[0].plot(temp_stats[4,:], temp_stats_minmax[0,:], label='min')
    axs[0].plot(temp_stats[4,:], temp_stats_minmax[1,:], label='max')
    axs[0].plot(temp_stats[4,:], temp_stats[0,:])
    axs[1].plot(temp_stats[4,:], temp_stats_minmax[2,:], label='min')
    axs[1].plot(temp_stats[4,:], temp_stats_minmax[3,:], label='max')
    axs[1].plot(temp_stats[4,:], temp_stats[1,:])
    axs[2].plot(temp_stats[4,:], temp_stats_minmax[4,:], label='min')
    axs[2].plot(temp_stats[4,:], temp_stats_minmax[5,:], label='max')
    axs[2].plot(temp_stats[4,:], temp_stats[2,:])
    axs[3].plot(temp_stats[4,:], temp_stats_minmax[6,:], label='min')
    axs[3].plot(temp_stats[4,:], temp_stats_minmax[7,:], label='max')
    #axs[3].plot(temp_stats[4,:], temp_stats[3,:])
    axs[3].plot(temp_stats[4,1:], temp_stats[3,1:])
    for i in range(4):
        #axs[i].set_yscale('log')
        axs[i].set_xlabel('Time')
        axs[i].set_ylabel('L2-norm')
        axs[i].grid(alpha=0.3, which='both', linestyle=':')
        axs[i].legend(loc="upper left")
    axs[0].set_title('$||c_a||_{L^2(Y)}$')
    axs[1].set_title('$||c_d||_{L^2(Y)}$')        
    axs[2].set_title('$||p_a||_{L^2(\Gamma)}$')
    axs[3].set_title('$||u||_{L^2(Y)}$')
    #plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.png')
    plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.png')
    #plt.show()
    plt.close()

    # plot the temporal statistics
    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten()
    axs[0].plot(variables[10,:], variables[0,:], label='min')
    axs[0].plot(variables[10,:], variables[1,:], label='max')
    axs[0].plot(variables[10,:], variables[2,:])
    axs[0].legend(loc="lower right")
    axs[1].plot(variables[10,:], variables[3,:], label='min')
    axs[1].plot(variables[10,:], variables[4,:], label='max')
    axs[1].plot(variables[10,:], variables[5,:])
    axs[1].legend(loc="upper left")
    axs[2].plot(variables[10,:], variables[6,:])
    axs[3].plot(variables[10,:], variables[7,:], label='min')
    axs[3].plot(variables[10,:], variables[8,:], label='max')
    axs[3].plot(variables[10,:], variables[9,:])
    axs[3].legend(loc="upper left")
    for i in range(4):
        #axs[i].set_yscale('log')
        axs[i].set_xlabel('Time')
        axs[i].grid(alpha=0.3, which='both', linestyle=':')
    axs[0].set_title('$||tr(sigma(u))||_{L^2(Y)}$')
    axs[1].set_title('$||E_c||_{L^2(Y)}$')
    axs[2].set_title('$||\partial_t div(u)||_{L^2(Y)}$')
    axs[3].set_title('$||div(u)||_{L^2(Y)}$')
    #plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'-oth-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.png')
    plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'-oth-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.png')
    #plt.show()
    plt.close()


def plot2D(coupling, meshnr, dt, T, k6, bigE, str_C1, bcs):

    temp_der_stats = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'_2D-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')
    temp_stats_minmax = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'_2D-minmax-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')
    temp_stats = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'_2D-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')
    variables = np.load('../../../new_results/temp/tempstats'+str(coupling)+str(bcs)+'_2D-oth-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.npy')


    # plot the temporal statistics
    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten() 
    axs[0].plot(temp_der_stats[4,:], temp_der_stats[0,:])
    axs[1].plot(temp_der_stats[4,:], temp_der_stats[1,:])
    axs[2].plot(temp_der_stats[4,:], temp_der_stats[2,:])
    #axs[3].plot(temp_der_stats[4,:], temp_der_stats[3,:])
    axs[3].plot(temp_der_stats[4,1:], temp_der_stats[3,1:])
    for i in range(4):
        #axs[i].set_yscale('log')
        axs[i].set_xlabel('Time')
        axs[i].set_ylabel('L2-norm')
        axs[i].grid(alpha=0.3, which='both', linestyle=':')
    axs[0].set_title('$||\partial_t c_a||_{L^2(Y)}$')
    axs[1].set_title('$||\partial_t c_d||_{L^2(Y)}$')
    axs[2].set_title('$||\partial_t p_a||_{L^2(\Gamma)}$')
    axs[3].set_title('$||\partial_t u||_{L^2(Y)}$')
    #plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'_2D-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+'_'+str(bigE)+'E.png')
    plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'_2D-der-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.png')
    #plt.show()
    plt.close()

    # plot the temporal statistics
    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten()
    axs[0].plot(temp_stats[4,:], temp_stats_minmax[0,:], label='min')
    axs[0].plot(temp_stats[4,:], temp_stats_minmax[1,:], label='max')
    axs[0].plot(temp_stats[4,:], temp_stats[0,:])
    axs[1].plot(temp_stats[4,:], temp_stats_minmax[2,:], label='min')
    axs[1].plot(temp_stats[4,:], temp_stats_minmax[3,:], label='max')
    axs[1].plot(temp_stats[4,:], temp_stats[1,:])
    axs[2].plot(temp_stats[4,:], temp_stats_minmax[4,:], label='min')
    axs[2].plot(temp_stats[4,:], temp_stats_minmax[5,:], label='max')
    axs[2].plot(temp_stats[4,:], temp_stats[2,:])
    axs[3].plot(temp_stats[4,:], temp_stats_minmax[6,:], label='min')
    axs[3].plot(temp_stats[4,:], temp_stats_minmax[7,:], label='max')
    #axs[3].plot(temp_stats[4,:], temp_stats[3,:])
    axs[3].plot(temp_stats[4,1:], temp_stats[3,1:])
    for i in range(4):
        #axs[i].set_yscale('log')
        axs[i].set_xlabel('Time')
        axs[i].set_ylabel('L2-norm')
        axs[i].grid(alpha=0.3, which='both', linestyle=':')
        axs[i].legend(loc="upper left")
    axs[0].set_title('$||c_a||_{L^2(Y)}$')
    axs[1].set_title('$||c_d||_{L^2(Y)}$')        
    axs[2].set_title('$||p_a||_{L^2(\Gamma)}$')
    axs[3].set_title('$||u||_{L^2(Y)}$')
    plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'_2D-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.png')
    #plt.show()
    plt.close()

    # plot the temporal statistics
    fig, axs = plt.subplots(2, 2, figsize=(12,12))
    axs = axs.flatten()
    axs[0].plot(variables[10,:], variables[0,:], label='min')
    axs[0].plot(variables[10,:], variables[1,:], label='max')
    axs[0].plot(variables[10,:], variables[2,:])
    axs[0].legend(loc="lower right")
    axs[1].plot(variables[10,:], variables[3,:], label='min')
    axs[1].plot(variables[10,:], variables[4,:], label='max')
    axs[1].plot(variables[10,:], variables[5,:])
    axs[1].legend(loc="upper left")
    axs[2].plot(variables[10,:], variables[6,:])
    axs[3].plot(variables[10,:], variables[7,:], label='min')
    axs[3].plot(variables[10,:], variables[8,:], label='max')
    axs[3].plot(variables[10,:], variables[9,:])
    axs[3].legend(loc="upper left")
    for i in range(4):
        #axs[i].set_yscale('log')
        axs[i].set_xlabel('Time')
        axs[i].grid(alpha=0.3, which='both', linestyle=':')
    axs[0].set_title('$||tr(sigma(u))||_{L^2(Y)}$')
    axs[1].set_title('$||E_c||_{L^2(Y)}$')
    axs[2].set_title('$||\partial_t div(u)||_{L^2(Y)}$')
    axs[3].set_title('$||div(u)||_{L^2(Y)}$')
    plt.savefig('../../../new_results/figures/tempstats'+str(coupling)+str(bcs)+'_2D-oth-m'+str(meshnr)+'_dt='+str(dt)+'_T='+str(T)+'_k6='+str(k6)+str_C1+'_'+str(bigE)+'E.png')
    #plt.show()
    plt.close()


def compute_EOCs(max_cellsizes, errors):
    """
    compute EOCs given error norms
    """
    EOCs = np.zeros((len(errors),len(errors[0])-1))
    for j in range(len(errors)):
        for i in range(len(errors[0])-1):
            EOCs[j,i] = np.log(errors[j][i+1]/errors[j][i]) / np.log(max_cellsizes[i+1]/max_cellsizes[i])

    return EOCs

"""
def build_nullspace(V, x):
    #Function to build null space for 3D elasticity
    

    # Create list of vectors for null space
    nullspace_basis = [x.copy() for i in range(6)]

    # Build translational null space basis
    V.sub(0).dofmap().set(nullspace_basis[0], 1.0);
    V.sub(1).dofmap().set(nullspace_basis[1], 1.0);
    V.sub(2).dofmap().set(nullspace_basis[2], 1.0);

    # Build rotational null space basis
    V.sub(0).set_x(nullspace_basis[3], -1.0, 1);
    V.sub(1).set_x(nullspace_basis[3],  1.0, 0);
    V.sub(0).set_x(nullspace_basis[4],  1.0, 2);
    V.sub(2).set_x(nullspace_basis[4], -1.0, 0);
    V.sub(2).set_x(nullspace_basis[5],  1.0, 1);
    V.sub(1).set_x(nullspace_basis[5], -1.0, 2);

    for x in nullspace_basis:
        x.apply("insert")

    # Create vector space basis and orthogonalize
    basis = VectorSpaceBasis(nullspace_basis)
    basis.orthonormalize()

    return basis


def create_meshes(M, name = 'unit-sphere'):
    #Function that creates an M number of refinedments of the 'name' meshes from gmsh files
    
    # Create mesh and define function space
    meshes = []
    for i in range(M):
        meshes.append(Mesh()) #define mesh

    for i in range(1,M+1):
        # write xdvmf file from msh file (only need to do this once)
        msh = meshio.read('mesh/'+name+'{}.msh'.format(i)) #read msh file
        meshio.write('mesh/'+name+'{}.xdmf'.format(i), meshio.Mesh(points=msh.points, cells={"tetra": msh.cells_dict["tetra"]})) #convert to xdmf file

        with XDMFFile('mesh/'+name+'{}.xdmf'.format(i)) as xdmf: 
            xdmf.read(meshes[i-1]) #write xdmf file to mesh
        File('mesh/'+name+'{}.pvd'.format(i)).write(meshes[i-1]) #save mesh to new pvd file
    
    return np.array(meshes)
"""