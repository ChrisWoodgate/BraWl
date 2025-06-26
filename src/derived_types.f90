!> @file    derived_types.f90
!>
!> @brief   Module containing derived_types.
!>
!> @details 
!>
!> @author  C. D. Woodgate
!> @author  H. J. Naguszewski
!> @author  L. B. Partay
!>
!> @date    2019-2025
module derived_types

  use kinds
  use constants
  
  implicit none

  private

  public :: run_params, metropolis_params, ns_params, tmmc_params, wl_params

  !> @brief   Derived type for parameters specifying general simulation
  !>          parameters which are common to all sampling methods.
  !>
  !> @details run_params is a derived type specifying general simulation
  !>          parameters which are common to all sampling methods
  !>
  !> @author  C. D. Woodgate
  !>
  !> @date    2020-2025
  !>
  !> @param mode Paramater defining the mode of operation
  !> @param n_1 Number of atoms in direction 1 (x, for orthorhombic systems)
  !> @param n_2 Number of atoms in direction 1 (x, for orthorhombic systems)
  !> @param n_3 Number of atoms in direction 1 (x, for orthorhombic systems)
  !> @param n_basis Number of atoms in basis (Note: not user-specified, determined based on lattice type)
  !> @param n_species Number of chemical species
  !> @param n_atoms Number of atoms in the cell (Note: not user-specified, determined based on size and lattice type)
  !> @param static_seed Whether to use a static PRNG seed (for testing)
  !> @param lattice String naming the lattice type. Currently 'fcc', 'bcc', and 'simple_cubic' are supported.
  !> @param lattice_parameter Lattice constant (Note: only used for i/o of xyz files. Does not affect calculations internally)
  !> @param lattice_vectors (Developmental.) Vectors specifying the underlying lattice
  !> @param basis_vectors (Developmental.) Set of basis vectors
  !> @param species_names Set of names of chemical species (elements)
  !> @param species_concentrations Concentrations of each chemical species (cannot specify this and species_numbers)
  !> @param species_numbers Number of atoms of each chemical species (cannot specify this and species_concentrations)
  !> @param interaction_file Name of file containing atom-atom effective pair interactions
  !> @param interaction_range Range (in number of coordination shells) of the EPIs
  !> @param wc_range Number of shells (including the on-site values) for calculation of the WC ASRO parameters
  !> @param full_energy Function pointer to full Hamiltoninan implmentation used at runtime
  !> @param nbr_energy Function pointer to local Hamiltoninan implmentation used at runtime
  !> @param rdm_site Function pointer to random lattice site implmentation used at runtime
  !> @param rdm_nbr Function pointer to random neighbour implmentation used at runtime
  !> @param mc_step Function pointer to Monte Carlo step implmentation used at runtime
  type run_params

    ! Paramater defining the mode of operation
    integer :: mode
    ! Number of unit cells in each direction and
    ! number of atoms in the basis
    integer :: n_1, n_2, n_3, n_basis
    ! Number of chemical species
    integer :: n_species
    ! Number of atoms
    integer :: n_atoms
    ! Logical flag for fixed or time-based random seed
    logical :: static_seed=.false.
    ! Lattice type - name, e.g. fcc, bcc, hcp, ...
    character(len=20) :: lattice
    ! Lattice parameter (for writing xyz file)
    real(real64) :: lattice_parameter
    ! Lattice vectors (for writing xyz file)
    real(real64), dimension(3,3) :: lattice_vectors
    ! Vector to second basis atom (for writing xyz file)
    real(real64), dimension(3) :: basis_vectors
    ! Names of the chemical species
    character(len=2), dimension(:), allocatable :: species_names
    ! Concentrations of the chemical species
    real(real64), dimension(:), allocatable :: species_concentrations
    ! Number of atoms of each chemical species
    integer(int64), dimension(:), allocatable :: species_numbers
    ! Atom-atom interchange interaction file name
    character(len=50) :: interaction_file
    ! Interaction range (number of coordination shells)
    integer :: interaction_range
    ! Number of coordination shells for Warren-Cowley parameters
    integer :: wc_range
    ! Hamiltonian to use
    procedure(hamiltonian), pointer, pass :: full_energy => null()
    ! neighbour energy function to use
    procedure(neighbour), pointer, pass :: nbr_energy => null()
    ! Random site function to use
    procedure(rand_site), pointer, pass :: rdm_site => null()
    ! Random neighbour function to use
    procedure(rand_neighbour), pointer, pass :: rdm_nbr => null()
    ! Monte Carlo step to call. (Neighbour swap or whole lattice swap)
    procedure(monte_carlo), pointer, pass :: mc_step => null()

  end type run_params

  !> @brief   Derived type for parameters defining a Metropolis-Hastings
  !>          run
  !>
  !> @details metropolis_params is a derived type specifying parameters
  !>          defining a Metropolis-Hastings run
  !>
  !> @author  C. D. Woodgate
  !>
  !> @date    2020-2025
  !>
  !> @param mode Paramater defining the mode of Metropolis operation
  !> @param n_mc_steps
  !> @param burn_in_start Whether to burn in at the initial temperature
  !> @param burn_in Whether to burn in at all temperatures
  !> @param n_burn_in_steps Number of trial Monte Carlo steps to use for burn-in
  !> @param n_sample_steps Number of trial Monte Carlo steps between collecting statistics
  !> @param calculate_energies Whether to output energy data. Defaults to True
  !> @param write_trajectory_energy Whether to output the energy as a function of step #. Default to False
  !> @param calculate_asro Whether to calculate ASRO parameters
  !> @param write_trajectory_asro Whether to output the ASRO as a function of step #. Default to False
  !> @param n_sample_steps_asro Number of trial Monte Carlo steps between sampling ASRO.
  !>                            Note: must be a multiple of n_sample_steps
  !> @param calculate_alro Whether to calculate ALRO parameters
  !> @param n_sample_steps_alro Number of trial Monte Carlo steps between sampling ALRO.
  !>                            Note: must be a multiple of n_sample_steps
  !> @param write_trajectory_xyz Whether to output the trajectory. Default to False---can generate BIG files.
  !> @param n_sample_steps_trajectory Number of trial Monte Carlo steps between sampling trajectory. 
  !>                                  Note: must be a multiple of n_sample_steps
  !> @param write_final_config_nc Whether to output the final config at each T in binary format
  !> @param write_final_config_xyz Whether to output the final config at each T in xyz format
  !> @param read_start_config_nc Whether to real the initial config from binary file
  !> @param start_config_file Name of file from which to read initial config
  !> @param beta Inverse temperature, internally in Ry
  !> @param T Temperature (or start temperature, if using simulated annealing)
  !> @param T_steps Number of temperature steps for simulated annealing
  !> @param delta_T Temperature step size for simulated annealing
  !> @param nbr_swap Whether to trial just neighbour swaps, or swaps across the whole lattice
  type metropolis_params

    ! Paramater defining the mode of Metropolis operation
    ! Current options are:
    !  - 'equilibrate'
    !    - This option runs the Metropolis algorithm at a given
    !      temperature and tracks the simulation energy and atomic
    !      order parameters as a function of number of trial moves.
    !      Optionally, you can dump a 'trajectory' showing how the
    !      alloy configuration evolves.
    !  - 'simulated_annealing'
    !    - This option uses simulated annealing to examine atomic
    !      short- and long-range order parameters, simulation energy,
    !      and specific heat capacity as a function of temperature for a
    !      given system. Optionally, you can dump the simulation
    !      configuration at the end of each temperature.
    !  - 'decorrelated_samples'
    !    - This option draws grid configurations a set number of steps
    !      apart (for producing indicative equilibrated atomic configs
    !      at a given temperature.
    character(len=20) :: mode
    ! Number of monte carlo steps (at each temperature if annealing)
    integer :: n_mc_steps
    ! Burn in at FIRST temperature if doing simulated annealing?
    logical :: burn_in_start
    ! Burn in at ALL temperatures if doing simulated annealing?
    logical :: burn_in
    ! Number of burn-in steps (at each temperature if annealing)
    integer :: n_burn_in_steps
    ! Number of Monte Carlo trial moves between drawing samples
    ! of ANY quantity
    integer :: n_sample_steps
    ! Do we want to total simulation energy and store it?
    logical :: calculate_energies
    ! Do we want to write the trajectory's energy to file
    logical :: write_trajectory_energy
    ! Do we want to calculate atomic short-range order parameters?
    logical :: calculate_asro
    ! Do we want to write the trajectory's energy to file
    logical :: write_trajectory_asro
    ! Number of monte carlo steps between drawing ASRO parameters
    ! (MUST be an integer multiple of n_sample_steps)
    integer :: n_sample_steps_asro
    ! Do we want to store atomic long-range order parameters
    logical :: calculate_alro
    ! Number of monte carlo steps between drawing ALRO parameters
    ! (MUST be an integer multiple of n_sample_steps)
    integer :: n_sample_steps_alro
    ! Do we want to write the trajectory as an xyz file
    ! (for visualisation)
    logical :: write_trajectory_xyz
    ! Number of monte carlo steps between drawing ALRO parameters
    ! (MUST be an integer multiple of n_sample_steps)
    integer :: n_sample_steps_trajectory
    ! Do we want to write the final config (at each temp if annealing)
    ! as an xyz file?
    logical :: write_final_config_xyz
    ! Do we want to write the final config (at each temp if annealing)
    ! as a NetCDF file? (For restart later.)
    logical :: write_final_config_nc
    ! Do we want to read a starting config from file?
    ! (Useful if restarting an earlier run...)
    logical :: read_start_config_nc
    ! Name of file from which to read if so
    character(len=144) :: start_config_file
    ! Inverse temperature
    real(real64) :: beta
    ! Temperature of simulation (or start temperature if annealing)
    ! Should come in units of Kelvin
    real(real64) :: T
    ! Number of temperature steps (if annealing)
    integer :: T_steps
    ! Temperature step size if annealing
    real(real64) :: delta_T
    ! Atom swap flag
    ! If True then only swap neighbouring pairs of atoms
    ! If False then permit swaps across the lattice
    logical :: nbr_swap

  end type metropolis_params

  !> @brief   Derived type for nested sampling parameters.
  !>
  !> @details ns_params is a structure consiting various data for nested sampling.
  !>
  !> @author  L. B. Partay
  !>
  !> @date    2024
  !>
  !> @param  n_walkers Number of nested sampling walkers
  !> @param  n_steps Number of steps required for walk generating a new configuration
  !> @param  n_iter Number of NS iterations before sampling is finished
  !> @param  traj_freq Frequency of writing configuration (every traj_freq-th NS iteration)
  !> @param  outfile_ener Output filename for energies
  !> @param  outfile_traj Output filename for configurations
  type ns_params

    ! Number of nested sampling walkers
    integer :: n_walkers
    ! Number of steps required for walk generating a new configuration
    integer :: n_steps
    ! Number of NS iterations before sampling is finished
    integer :: n_iter
    ! Frequency of writing configuration (every traj_freq-th NS iteration)
    integer :: traj_freq
    ! Output filenames
    character(len=100) :: outfile_ener, outfile_traj  

  end type ns_params

  !> @brief   Derived type for TMMC parameters.
  !>
  !> @details tmmc_params is a structure consiting various data for TMMC.
  !>
  !> @author  H. J. Naguszewski
  !>
  !> @date    2024
  !>
  !> @param  mc_sweeps Number of sweeps (each sweep is n_atoms mc steps)
  !> @param  bins Number of bins across energy range
  !> @param  num_windows Number of energy windows
  !> @param  bin_overlap Number of bins in the overlap region
  !> @param  weight_update Number of bias weight updates
  !> @param  energy_min Energy range minimum
  !> @param  energy_max Energy range maximum
  !> @param  T Temperature to use for dynamics
  type tmmc_params

    ! Number of sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
    ! Number of bins across energy range
    integer :: bins
    ! Number of energy windows
    integer :: num_windows
    ! Number of bins in the overlap region
    real :: bin_overlap
    ! Number of bias weight updates
    integer :: weight_update
    ! Energy range minimum
    real :: energy_min
    ! Energy range maximum
    real :: energy_max
    ! Temperature of dynamics
    real :: T

  end type tmmc_params

  !> @brief   Derived type for Wang-Landau sampling parameters.
  !>
  !> @details wl_params is a derived type consiting various data for
  !>          Wang-Landau sampling.
  !>
  !> @author  H. J. Naguszewski
  !>
  !> @date    2024
  !>
  !> @param  mc_sweeps Number of sweeps (each sweep is n_atoms mc steps)
  !> @param  bins Number of bins across energy range
  !> @param  num_windows Number of energy windows
  !> @param  bin_overlap Percentage overlap between windows
  !> @param  tolerance Tolerance for Wang-Landau
  !> @param  wl_f Flatness for Wang-Landau histogram
  !> @param  energy_min Energy range minimum
  !> @param  energy_max Energy range maximum
  !> @param  T Temperature to use for dynamics
  !> @param  radial_samples Number of radial density samples to draw per bin
  !> @param  performance Performance analysis mode 
  !>         0 - non-uniform windows | dynamic window sizes | replica exchange
  !>         1 - non-uniform windows | dynamic window sizes |        x
  !>         2 - non-uniform windows |          x           | replica exchange
  !>         3 - non-uniform windows |          x           |        x
  !>         4 -         x           |          x           | replica exchange
  !>         5 -         x           |          x           |        x
  type wl_params

    ! Number of  sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
    ! Number of bins across energy range
    integer :: bins
    ! Number of energy windows
    integer :: num_windows
    ! Percentage overlap between windows
    real :: bin_overlap
    ! Tolerance for wang landau
    real :: tolerance
    ! Flatness for wang landau histogram
    real :: flatness
    ! Wang Landau density of states histogram tuning parameter
    real :: wl_f
    ! Energy range minimum
    real :: energy_min
    ! Energy range maximum
    real :: energy_max
    ! Radial density samples per bin
    integer :: radial_samples
    ! Performance analysis
    integer :: performance

  end type wl_params

  !--------------------------------------------------------------------!
  ! Interfaces for various Hamiltonian and dynamics implementations    !
  ! chosen at runtime.                                                 !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  interface

    !> @brief   Function for evaluating the energy of the simulation
    !>
    !> @details This is set to 'point' at the relevant implementation at
    !>          runtime.
    !>
    !> @author  C. D. Woodgate
    !>
    !> @date    2020-2025
    !>
    !> @param setup Derived type specifying the simulation
    !> @param config The configuration for which to evaluate the energy
    !>
    !> @return The simulation energy
    function hamiltonian(setup, config)

      use kinds, only : real64
      use constants, only : array_int

      import :: run_params

      integer(array_int), dimension(:,:,:,:), intent(in) :: config
      real(real64) :: hamiltonian
      class(run_params), intent(in) :: setup

    end function

    !> @brief   Function for getting a random neighbour of a lattice
    !>          site
    !>
    !> @author  C. D. Woodgate
    !>
    !> @date    2020-2025
    !>
    !> @param setup Derived type specifying the simulation
    !> @param config The configuration for which to evaluate the energy
    !> @param site_b Site basis index
    !> @param site_1 Site index 1
    !> @param site_2 Site index 2
    !> @param site_3 Site index 3
    !>
    !> @return An integer indexing the neighbour
    function neighbour(setup, config, site_b, site_i, site_j, site_k)

      use kinds, only : real64
      use constants, only : array_int

      import :: run_params

      integer(array_int), dimension(:,:,:,:), intent(in) :: config
      real(real64) :: neighbour
      class(run_params), intent(in) :: setup
      integer, intent(in) :: site_b, site_i, site_j, site_k

    end function

    !> @brief   Function for getting a random lattice site
    !>
    !> @author  C. D. Woodgate
    !>
    !> @date    2020-2025
    !>
    !> @param setup Derived type specifying the simulation
    !>
    !> @return Array of four ints specifying the site
    function rand_site(setup)

      use kinds, only : real64
      use constants, only : array_int

      import :: run_params

      class(run_params), intent(in) :: setup
      integer, dimension(4) :: rand_site

    end function

    !> @brief   Function for getting a random neighbour of a lattice
    !>          site
    !>
    !> @author  C. D. Woodgate
    !>
    !> @date    2020-2025
    !>
    !> @param setup Derived type specifying the simulation
    !> @param site Array of four ints specifying the site
    !>
    !> @return Array of four ints specifying the neighbour site
    function rand_neighbour(setup, site)

      use kinds, only : real64
      use constants, only : array_int

      import :: run_params

      class(run_params), intent(in) :: setup
      integer, dimension(4), intent(in) :: site
      integer, dimension(4) :: rand_neighbour

    end function

    !> @brief   Function performing a trial Monte Carlo move on a config
    !>
    !> @details This is set to 'point' at the relevant implementation at
    !>          runtime.
    !>
    !> @author  C. D. Woodgate
    !>
    !> @date    2020-2025
    !>
    !> @param setup Derived type specifying the simulation
    !> @param config The current configuration
    !> @param config The current configuration
    !>
    !> @return 1 if the trial move is accepted, 0 otherwise
    function monte_carlo(setup, config, beta) result(accept)

      use kinds, only : real64
      use constants, only : array_int

      import :: run_params

      integer(array_int), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      integer :: accept
      real(real64), intent(in) :: beta

    end function monte_carlo

  end interface

end module derived_types
