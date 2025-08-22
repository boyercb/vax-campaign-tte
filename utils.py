import random 
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from lifelines import CoxPHFitter
from sklearn.linear_model import LogisticRegression

def generate_random_graph(num_nodes, edges_per_node):
    """Generates a random graph with a specified number of nodes and expected edges."""
    probability = edges_per_node / (num_nodes - 1)
    G = nx.erdos_renyi_graph(num_nodes, probability)
    return G

def assign_initial_states(G, infection_prob, recovery_prob, death_prob, testing_prob, vaccine_efficacy):
    """Assigns initial states to nodes in the graph."""
    for node in G.nodes():
        G.nodes[node]['S'] = True   # Susceptible
        G.nodes[node]['I'] = False  # Infectious
        G.nodes[node]['R'] = False  # Recovered
        G.nodes[node]['D'] = False  # Dead
        G.nodes[node]['T'] = False  # Tested
        G.nodes[node]['U'] = random.random()  # Unmeasured confounder
        G.nodes[node]['X'] = random.random()  # Measured confounder
        G.nodes[node]['time'] = 0   # Follow-up time (for censoring)
        G.nodes[node]['infection_time'] = 0  # Time of infection
        G.nodes[node]['test_time'] = 0       # Time of testing
        G.nodes[node]['death_time'] = 0      # Time of death
        G.nodes[node]['vaccinated'] = False  # Not vaccinated initially
        G.nodes[node]['infection_prob'] = infection_prob  # Base infection probability
        G.nodes[node]['recovery_prob'] = recovery_prob    # Probability of recovery
        G.nodes[node]['death_prob'] = death_prob          # Probability of death
        G.nodes[node]['testing_prob'] = testing_prob      # Probability of testing
        G.nodes[node]['vaccine_efficacy'] = vaccine_efficacy  # Efficacy of vaccine
        
        # Calculate individual infection probabilities based on confounders
        G.nodes[node]['infection_prob'] *= np.exp(G.nodes[node]['X'] - 0.5) 
        G.nodes[node]['infection_prob_U'] = np.exp(G.nodes[node]['U'] - 0.5) * G.nodes[node]['infection_prob']

        # Calculate individual testing probabilities based on confounders
        G.nodes[node]['testing_prob'] *= np.exp(G.nodes[node]['X'] - 0.5) 
        G.nodes[node]['testing_prob_U'] = np.exp(G.nodes[node]['U'] - 0.5) * G.nodes[node]['testing_prob']

def assign_initial_infected(G, num_initial_infected):
    """Assigns initial infected nodes in the graph."""
    initial_infected = random.sample(G.nodes(), num_initial_infected)
    for node in initial_infected:
        G.nodes[node]['I'] = True
        G.nodes[node]['S'] = False
        G.nodes[node]['infection_time'] = 0  # Initially infected at time 0

def vaccination_campaign(G, base_vaccination_rate, step):
    """
    Applies vaccination campaign at each time step to unvaccinated individuals.
    Vaccination probability varies by individual characteristics X and U.
    
    Parameters:
    - G: Graph with node attributes
    - base_vaccination_rate: Base proportion to vaccinate at this step
    - step: Current time step
    """
    for node in G.nodes():
        # Only vaccinate those who are not already vaccinated and are alive
        if not G.nodes[node]['vaccinated'] and not G.nodes[node]['D']:
            # Individual vaccination probability based on confounders
            individual_vax_prob = base_vaccination_rate * np.exp(G.nodes[node]['U'] + G.nodes[node]['X'] - 1)
            
            if random.random() < individual_vax_prob:
                G.nodes[node]['vaccinated'] = True
                
                # Adjust infection probabilities based on vaccine efficacy
                G.nodes[node]['infection_prob'] *= (1 - G.nodes[node]['vaccine_efficacy'])
                G.nodes[node]['infection_prob_U'] *= (1 - G.nodes[node]['vaccine_efficacy'])

def simulate_outbreak(G, steps, plot=True, print_progress=True, vaccination_schedule=None):
    """
    Simulates disease outbreak in the graph over a number of steps.
    
    Parameters:
    - G: Graph with node attributes
    - steps: Number of simulation steps
    - plot: Whether to plot epidemic curves
    - print_progress: Whether to print progress updates
    - vaccination_schedule: List of vaccination rates for each step (optional)
    """
    # Initialize tracking arrays
    infected_over_time = []
    recovered_over_time = []
    dead_over_time = []
    tested_over_time = []
    vaccinated_over_time = []
    tested_over_time_vaccinated = []
    tested_over_time_unvaccinated = []
    
    for step in range(steps):
        
        # Apply vaccination campaign if schedule is provided
        if vaccination_schedule and step < len(vaccination_schedule):
            vaccination_campaign(G, vaccination_schedule[step], step)
        
        for node in G.nodes():
            # If infected
            if G.nodes[node]['I']:
                # Check for death first
                if random.random() < G.nodes[node]['death_prob']:
                    G.nodes[node]['D'] = True
                    G.nodes[node]['I'] = False
                    G.nodes[node]['death_time'] = step
                # Otherwise check for recovery
                elif random.random() < G.nodes[node]['recovery_prob']:
                    G.nodes[node]['R'] = True
                    G.nodes[node]['I'] = False
                else:
                    # Spread infection to neighbors if they are susceptible
                    for neighbor in G.neighbors(node):
                        if G.nodes[neighbor]['S'] and not G.nodes[neighbor]['D']:
                            if step < 20:  # Only measured confounder in the first 20 steps
                                if random.random() < G.nodes[neighbor]['infection_prob']:
                                    G.nodes[neighbor]['I'] = True
                                    G.nodes[neighbor]['S'] = False
                                    G.nodes[neighbor]['infection_time'] = step
                                    if random.random() < G.nodes[neighbor]['testing_prob']:
                                        G.nodes[neighbor]['T'] = True
                                        G.nodes[neighbor]['test_time'] = step
                            else:  # After 20 steps, activate unmeasured confounder
                                if random.random() < G.nodes[neighbor]['infection_prob_U']:
                                    G.nodes[neighbor]['I'] = True
                                    G.nodes[neighbor]['S'] = False
                                    G.nodes[neighbor]['infection_time'] = step
                                    if random.random() < G.nodes[neighbor]['testing_prob_U']:
                                        G.nodes[neighbor]['T'] = True
                                        G.nodes[neighbor]['test_time'] = step

            # Update follow-up time for all living individuals
            if not G.nodes[node]['D']:
                G.nodes[node]['time'] = step

        # Track epidemic curve by calculating the number of nodes in each state
        num_infected = sum(1 for n in G.nodes() if G.nodes[n]['I'])
        num_recovered = sum(1 for n in G.nodes() if G.nodes[n]['R'])
        num_dead = sum(1 for n in G.nodes() if G.nodes[n]['D'])
        num_tested = sum(1 for n in G.nodes() if G.nodes[n]['T'])
        num_vaccinated = sum(1 for n in G.nodes() if G.nodes[n]['vaccinated'])
        num_tested_vaccinated = sum(1 for n in G.nodes() if G.nodes[n]['T'] and G.nodes[n]['vaccinated'])
        num_tested_unvaccinated = sum(1 for n in G.nodes() if G.nodes[n]['T'] and not G.nodes[n]['vaccinated'])
        
        infected_over_time.append(num_infected)
        recovered_over_time.append(num_recovered)
        dead_over_time.append(num_dead)
        tested_over_time.append(num_tested)
        vaccinated_over_time.append(num_vaccinated)
        tested_over_time_vaccinated.append(num_tested_vaccinated)
        tested_over_time_unvaccinated.append(num_tested_unvaccinated)

        # Print progress
        if step % 10 == 0 and print_progress:
            print(f"Step {step}: Infected: {num_infected}, Recovered: {num_recovered}, Dead: {num_dead}, Tested: {num_tested}, Vaccinated: {num_vaccinated}")

    # Plot curves at the end
    if plot:
        # Plot epidemic curve
        plt.figure(figsize=(15, 10))
        
        plt.subplot(2, 2, 1)
        plt.plot(range(steps), infected_over_time, label='Infected', color='red')
        plt.plot(range(steps), recovered_over_time, label='Recovered', color='green')
        plt.plot(range(steps), dead_over_time, label='Dead', color='black')
        plt.xlabel('Time Steps')
        plt.ylabel('Number of Individuals')
        plt.title('Epidemic Curve')
        plt.legend()
        
        # Plot testing curve
        plt.subplot(2, 2, 2)
        plt.plot(range(steps), tested_over_time, label='Total Tested', color='blue')
        plt.xlabel('Time Steps')
        plt.ylabel('Number of Tested Individuals')
        plt.title('Testing Curve')
        plt.legend()
        
        # Plot vaccination curve
        plt.subplot(2, 2, 3)
        plt.plot(range(steps), vaccinated_over_time, label='Vaccinated', color='purple')
        plt.xlabel('Time Steps')
        plt.ylabel('Number of Vaccinated Individuals')
        plt.title('Vaccination Curve')
        plt.legend()
        
        # Plot testing curve by vaccination status
        plt.subplot(2, 2, 4)
        plt.plot(range(steps), tested_over_time_vaccinated, label='Tested Vaccinated', color='purple')
        plt.plot(range(steps), tested_over_time_unvaccinated, label='Tested Unvaccinated', color='orange')
        plt.xlabel('Time Steps')
        plt.ylabel('Number of Tested Individuals')
        plt.title('Testing by Vaccination Status')
        plt.legend()
        
        plt.tight_layout()
        plt.show()

    return G

def estimate_models(G, outcome='tested', print_results=True):
    """
    Estimates models based on the graph's node attributes.
    
    Parameters:
    - G: Graph with node attributes
    - outcome: Outcome to analyze ('tested', 'infected', or 'death')
    - print_results: Whether to print results
    
    Returns:
    - hr_basic: Hazard ratio without unmeasured confounder
    - hr_with_u: Hazard ratio with unmeasured confounder
    """
    # Create a DataFrame from the graph's node attributes
    data = pd.DataFrame.from_dict(dict(G.nodes(data=True)), orient='index')
    
    # Create event and time variables based on outcome
    if outcome == 'tested':
        data['event'] = data['T'].astype(int)
        # For tested individuals, use test time; for untested, use follow-up time
        data['event_time'] = np.where(data['T'], data['test_time'], data['time'])
    elif outcome == 'infected':
        # Create infection event (ever infected = not susceptible and not initially infected)
        data['event'] = (~data['S'] & (data['I'] | data['R'] | data['D'])).astype(int)
        # For infected individuals, use infection time; for non-infected, use follow-up time
        data['event_time'] = np.where(data['event'], data['infection_time'], data['time'])
    elif outcome == 'death':
        data['event'] = data['D'].astype(int)
        # For dead individuals, use death time; for alive, use follow-up time
        data['event_time'] = np.where(data['D'], data['death_time'], data['time'])
    else:
        raise ValueError("Outcome must be 'tested', 'infected', or 'death'")

    # Fit Cox proportional hazards model
    df = data[['event_time', 'event', 'vaccinated', 'X']].copy()
    model = CoxPHFitter()
    model.fit(df, duration_col='event_time', event_col='event')
    
    # Fit model with unmeasured confounder
    df_with_u = data[['event_time', 'event', 'vaccinated', 'X', 'U']].copy()
    model_with_u = CoxPHFitter()
    model_with_u.fit(df_with_u, duration_col='event_time', event_col='event')

    # Extract hazard ratios for vaccination status
    hr = np.exp(model.params_)
    hr_with_u = np.exp(model_with_u.params_)
    
    if print_results:
        print(f"Hazard Ratios for {outcome} (without unmeasured confounder):\n", hr)
        print(f"Hazard Ratios for {outcome} (with unmeasured confounder):\n", hr_with_u)

    # Return hazard ratio estimates for vaccination status
    return hr['vaccinated'], hr_with_u['vaccinated']

def run_simulation(num_nodes=1000, edges_per_node=5, infection_prob=0.02, 
                  recovery_prob=0.1, death_prob=0.01, testing_prob=0.3, 
                  vaccine_efficacy=0.8, num_initial_infected=5, steps=100,
                  vaccination_schedule=None, plot=True, print_progress=True,
                  outcome='tested'):
    """
    Convenience function to run a complete simulation.
    
    Parameters:
    - num_nodes: Number of nodes in the graph
    - edges_per_node: Expected number of edges per node
    - infection_prob: Base probability of infection transmission
    - recovery_prob: Probability of recovery per time step
    - death_prob: Probability of death per time step
    - testing_prob: Base probability of being tested
    - vaccine_efficacy: Effectiveness of vaccine (0-1)
    - num_initial_infected: Number of initially infected individuals
    - steps: Number of simulation steps
    - vaccination_schedule: List of vaccination rates per step
    - plot: Whether to generate plots
    - print_progress: Whether to print progress updates
    - outcome: Outcome to analyze ('tested', 'infected', or 'death')
    
    Returns:
    - G: Final graph state
    - results: Dictionary with hazard ratio estimates
    """
    
    # Create graph and initialize
    G = generate_random_graph(num_nodes, edges_per_node)
    assign_initial_states(G, infection_prob, recovery_prob, death_prob, 
                         testing_prob, vaccine_efficacy)
    assign_initial_infected(G, num_initial_infected)
    
    # Run simulation
    G = simulate_outbreak(G, steps, plot=plot, print_progress=print_progress, 
                         vaccination_schedule=vaccination_schedule)
    
    # Estimate models
    hr_basic, hr_with_u = estimate_models(G, outcome=outcome, print_results=print_progress)
    
    results = {
        'hr_basic': hr_basic,
        'hr_with_unmeasured_confounder': hr_with_u,
        'outcome': outcome
    }
    
    return G, results
