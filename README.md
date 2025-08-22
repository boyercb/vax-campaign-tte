# Vaccination Campaign Target Trial Emulation

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)

## Overview

This repository contains code for simulating epidemic outbreaks with vaccination campaigns using network-based models. The simulation framework supports time-to-event analysis of vaccination effectiveness across multiple outcomes (infection, testing, death) while accounting for measured and unmeasured confounding.

**Key Features:**
- Network-based epidemic modeling using social graphs
- Time-varying vaccination campaigns with individual-level heterogeneity
- Multiple outcome analysis (infection, testing positive, death)
- Proper handling of event timing for survival analysis
- Monte Carlo simulation for uncertainty quantification
- Parameter sensitivity analysis

## Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager

### Setup

1. **Clone the repository:**
   ```bash
   git clone https://github.com/[your-username]/vax-campaign-tte.git
   cd vax-campaign-tte
   ```

2. **Create a virtual environment (recommended):**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install required packages:**
   ```bash
   pip install -r requirements.txt
   ```

## Required Dependencies

The simulation requires the following Python packages:
- `numpy` - Numerical computing
- `pandas` - Data manipulation and analysis
- `matplotlib` - Plotting and visualization
- `networkx` - Network graph creation and analysis
- `lifelines` - Survival analysis (Cox proportional hazards models)
- `scikit-learn` - Machine learning utilities
- `seaborn` - Statistical data visualization
- `tqdm` - Progress bars
- `jupyter` - Jupyter notebook support

## Quick Start

### Basic Simulation

Run a basic epidemic simulation with vaccination campaign:

```python
from utils import run_simulation

# Define vaccination schedule (5% daily for first 50 days)
vaccination_schedule = [0.05] * 50 + [0.0] * 50

# Run simulation
G, results = run_simulation(
    num_nodes=1000,
    vaccine_efficacy=0.8,
    vaccination_schedule=vaccination_schedule,
    outcome='tested'
)

print(f"Vaccine effectiveness: {(1 - results['hr_basic']) * 100:.1f}%")
```

### Interactive Analysis

For comprehensive analysis, open the Jupyter notebook:

```bash
jupyter notebook simulation_notebook.ipynb
```

The notebook includes:
- Complete simulation examples
- Network visualization
- Vaccination strategy comparisons
- Parameter sensitivity analysis
- Monte Carlo uncertainty quantification

## Repository Structure

```
vax-campaign-tte/
├── utils.py                    # Core simulation functions
├── simulation_notebook.ipynb   # Interactive analysis notebook
├── README.md                   # This file
├── requirements.txt            # Python dependencies
├── LICENSE                     # MIT license
└── examples/                   # Additional example scripts (optional)
```

## Core Functions

### `run_simulation()`
Main function to execute a complete epidemic simulation with vaccination campaign.

**Key Parameters:**
- `num_nodes`: Population size (default: 1000)
- `vaccine_efficacy`: Vaccine effectiveness 0-1 (default: 0.8)
- `vaccination_schedule`: List of daily vaccination rates
- `outcome`: Analysis outcome ('tested', 'infected', 'death')

### `estimate_models()`
Estimates Cox proportional hazards models for vaccine effectiveness.

**Returns:**
- Hazard ratios with and without unmeasured confounding
- Proper event timing for each outcome type

### Network Functions
- `generate_random_graph()`: Creates Erdős–Rényi random graphs
- `assign_initial_states()`: Initializes node attributes and confounders
- `simulate_outbreak()`: Runs the epidemic simulation with vaccination

## Simulation Model

### Network Structure
- **Nodes**: Represent individuals in the population
- **Edges**: Represent social contacts through which disease spreads
- **Graph Type**: Erdős–Rényi random graph (configurable)

### Disease States
- **S**: Susceptible to infection
- **I**: Infectious (can transmit to neighbors)
- **R**: Recovered (immune)
- **D**: Dead (removed from transmission)
- **T**: Tested positive

### Confounding Structure
- **X**: Measured confounder (affects infection, testing, vaccination)
- **U**: Unmeasured confounder (affects infection, testing, vaccination)
- **Time-varying**: Unmeasured confounding activates after step 20

### Vaccination Campaign
- **Individual-level heterogeneity**: Vaccination probability varies by X and U
- **Time-varying rates**: Different vaccination rates per time step
- **Efficacy**: Reduces infection probability by specified amount

## Example Use Cases

### 1. Compare Vaccination Strategies
```python
schedules = {
    'Early': [0.05] * 50 + [0.0] * 50,
    'Late': [0.0] * 50 + [0.05] * 50,
    'Constant': [0.025] * 100
}

for name, schedule in schedules.items():
    G, results = run_simulation(vaccination_schedule=schedule, plot=False)
    print(f"{name}: VE = {(1-results['hr_basic'])*100:.1f}%")
```

### 2. Analyze Multiple Outcomes
```python
outcomes = ['infected', 'tested', 'death']
for outcome in outcomes:
    hr_basic, hr_confounded = estimate_models(G, outcome=outcome)
    print(f"{outcome.capitalize()}: HR = {hr_basic:.3f}")
```

### 3. Parameter Sensitivity
```python
efficacy_values = [0.4, 0.6, 0.8, 0.9]
for ve in efficacy_values:
    G, results = run_simulation(vaccine_efficacy=ve, plot=False)
    print(f"VE {ve*100}%: HR = {results['hr_basic']:.3f}")
```

## Output and Results

### Visualizations
- **Epidemic curves**: Disease progression over time
- **Vaccination curves**: Cumulative vaccination coverage
- **Network plots**: Final states and vaccination status
- **Survival curves**: Time-to-event analysis

### Statistical Output
- **Hazard ratios**: Vaccine effectiveness estimates
- **Confidence intervals**: From Monte Carlo simulation
- **Model comparisons**: With/without unmeasured confounding

<!-- ## Academic Citation

If you use this code in academic work, please cite:

```bibtex
@misc{vax_campaign_tte_2025,
  title={Vaccination Campaign Time-to-Event Analysis: Network-based Epidemic Simulation},
  author={[Your Name]},
  year={2025},
  publisher={GitHub},
  url={https://github.com/[your-username]/vax-campaign-tte}
}
``` -->

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests.

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/new-feature`)
3. Commit your changes (`git commit -am 'Add new feature'`)
4. Push to the branch (`git push origin feature/new-feature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or collaboration inquiries, please contact [boyerc5@ccf.org] or open an issue on GitHub.

## Acknowledgments

- Built using NetworkX for graph operations
- Survival analysis implemented with lifelines package
- Visualization powered by matplotlib and seaborn
