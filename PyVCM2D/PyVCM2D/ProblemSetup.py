"""
Module with functionalities to read simulation setup
"""
import yaml

class ReadYaml:

    def __init__(self, case = None):

        print('\n### ReadYaml class initialized ###')

        # Initialize attributes
        self.nx = None
        self.ny = None
        self.lx = None
        self.ly = None
        self.T_melt = None
        self.L_pcm = None
        self.cp_pcm = None
        self.k_pcm = None
        self.rho_pcm = None
        self.cp_bat = None
        self.k_bat = None
        self.rho_bat = None
        self.Q_bc = None
        self.T_bc = None
        self.T_initial = None
        self.total_time = None
        self.resmax = None
        self.kappa = None
        self.time_step = None

        # Read file
        if case is not None:
            self.case = case
            self.read_simulation_info()
        else:
            return

        # CV dimensions
        self.dx = self.lx/self.nx
        self.dy = self.ly/self.ny

        # Total number of CVs
        self.N = self.nx * self.ny

    def read_simulation_info(self):
        # Read file
        with open(self.case, 'r') as file:
            data = yaml.safe_load(file)

        # Get mesh info
        mesh_elements = data['Elements']
        self.nx = mesh_elements['x']
        self.ny = mesh_elements['y']
        mesh_length = data['Length']
        self.lx = mesh_length['x']
        self.ly = mesh_length['y']

        # Get physical properties
        battery = data['Battery']
        self.rho_bat = battery['rho']
        self.k_bat = battery['k']
        self.cp_bat = battery['cp']
        pcm = data['PCM']
        self.rho_pcm = pcm['rho']
        self.k_pcm = pcm['k']
        self.cp_pcm = pcm['cp']
        self.L_pcm = pcm['L']
        self.T_melt = pcm['T_melt']

        # Get initial and boundary conditions
        ic = data['Initial']
        self.T_initial = ic['Temp']
        bc = data['Boundary']
        self.T_bc = bc['Temp']
        self.Q_bc = bc['Flux']

        # Get simulation parameters
        self.time_step = data['Time_Step']
        self.kappa = data['Kappa']
        self.resmax = data['ResMax']
        self.total_time = data['Total_Time']
