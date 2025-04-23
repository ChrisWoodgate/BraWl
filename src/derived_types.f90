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

  public :: run_params, metropolis_params, ns_params, tmmc_params, &
            wl_params, es_params

  !--------------------------------------------------------------------!
  ! Type storing parameters defining simulation (used at runtime)      !
  !                                                                    !
  ! C. D. Woodgate,  Warwick                                      2023 !
  !--------------------------------------------------------------------!
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

  !--------------------------------------------------------------------!
  ! Type storing parameters defining Metropolis simulation             !
  ! (used at runtime)                                                  !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
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
    ! Burn in if doing simulated annealing?
    logical :: burn_in
    ! Number of burn-in steps (at each temperature if annealing)
    integer :: burn_in_steps
    ! Number of monte carlo steps (at each temperature if annealing)
    integer :: mc_steps
    ! Number of monte carlo steps between drawing energies
    integer :: sample_steps
    ! Number of monte carlo steps between drawing radial densities
    ! (MUST be a multiple of sample_steps)
    integer :: radial_sample_steps
    ! Do we want to store atomic short-range order parameters
    logical :: asro
    ! Do we want to store atomic long-range order parameters
    logical :: alro
    ! Do we want to store grids
    logical :: dump_grids
    ! Inverse temperature
    real(real64) :: beta
    ! Temperature of simulation (or start temperature if annealing)
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
  !> @param  n_walkers Number of nested sampling walkers
  !> @param  n_steps Number of steps required for walk generating a new configuration
  !> @param  n_iter Number of NS iterations before sampling is finished                   
  !> @param  traj_freq Frequency of writing configuration (every traj_freq-th NS iteration)
  !> @param  outfile_ener Output filename for energies
  !> @param  outfile_traj Output filename for configurations
  !> 
  !> @author  L. B. Partay
  !> @date    2024  
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

  !--------------------------------------------------------------------!
  ! Type storing parameters defining tmmc simulation (used at runtime) !
  !                                                                    !
  ! H. Naguszewski,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
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

  !--------------------------------------------------------------------!
  ! Type storing parameters defining wang landau simulation            !
  !  (used at runtime)                                                 !
  !                                                                    !
  ! H. Naguszewski,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
  type wl_params
    ! Number of  sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
    ! Number of bins across energy range
    integer :: bins
    ! Number of energy windows
    integer :: num_windows
    ! Number of bins in the overlap region
    integer :: bin_overlap
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

  end type wl_params

  !--------------------------------------------------------------------!
  ! Type storing parameters defining energy spectrum determining       !
  !  (used at runtime)                                                 !
  !                                                                    !
  ! H. Naguszewski,  Warwick                                      2024 !
  !--------------------------------------------------------------------!
  type es_params
    ! Number of mc sweeps (each sweep is n_atoms mc steps)
    integer :: mc_sweeps
  end type es_params

  !--------------------------------------------------------------------!
  ! Interfaces for various Hamiltonian and dynamics implementations    !
  ! chosen at runtime.                                                 !
  !                                                                    !
  ! C. D. Woodgate,  Bristol                                      2025 !
  !--------------------------------------------------------------------!
  interface

    ! Simulation energy
    function hamiltonian(setup, config)
      use kinds, only : real64
      use constants, only : array_int
      import :: run_params
      !integer(array_int), allocatable, dimension(:,:,:,:), intent(in) :: config
      integer(array_int), dimension(:,:,:,:), intent(in) :: config
      real(real64) :: hamiltonian
      class(run_params), intent(in) :: setup
    end function

    ! Neighbour
    function neighbour(setup, config, site_b, site_i, site_j, site_k)
      use kinds, only : real64
      use constants, only : array_int
      import :: run_params
      !integer(array_int), allocatable, dimension(:,:,:,:), intent(in) :: config
      integer(array_int), dimension(:,:,:,:), intent(in) :: config
      real(real64) :: neighbour
      class(run_params), intent(in) :: setup
      integer, intent(in) :: site_b, site_i, site_j, site_k
    end function

    ! Random site on the lattice
    function rand_site(setup)
      use kinds, only : real64
      use constants, only : array_int
      import :: run_params
      class(run_params), intent(in) :: setup
      integer, dimension(4) :: rand_site
    end function

    ! Random neighbour of that site
    function rand_neighbour(setup, site)
      use kinds, only : real64
      use constants, only : array_int
      import :: run_params
      class(run_params), intent(in) :: setup
      integer, dimension(4), intent(in) :: site
      integer, dimension(4) :: rand_neighbour
    end function

    ! Type of Monte Carlo step
    function monte_carlo(setup, config, beta) result(accept)
      use kinds, only : real64
      use constants, only : array_int
      import :: run_params
      !integer(array_int), allocatable, dimension(:,:,:,:) :: config
      integer(array_int), dimension(:,:,:,:) :: config
      class(run_params), intent(in) :: setup
      integer :: accept
      real(real64), intent(in) :: beta
    end function monte_carlo

  end interface

end module derived_types
