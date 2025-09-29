
# Multi-Domain Optimization Framework for ISAC

This repository contains the MATLAB implementation for the paper:

**R. Liu, M. Li, M. Zafari, B. Ottersten, and A. L. Swindlehurst, "Multi-domain optimization framework for ISAC: From electromagnetic shaping to network cooperation," *IEEE Wireless Communications*, under revision.**

## 📖 Abstract

Integrated sensing and communication (ISAC) has emerged as a key feature for sixth-generation (6G) networks, providing an opportunity to meet the dual demands of communication and sensing. Existing ISAC research primarily focuses on baseband optimization at individual access points, with limited attention to the roles of electromagnetic (EM) shaping and network-wide coordination. The intricate interdependencies between these domains remain insufficiently explored, leaving their full potential for enhancing ISAC performance largely untapped. To bridge this gap, we consider multi-domain ISAC optimization integrating EM shaping, baseband processing, and network cooperation strategies that facilitate efficient resource management and system-level design. We analyze the fundamental trade-offs between these domains and offer insights into domain-specific and cross-domain strategies contributing to ISAC performance and efficiency. We then conduct a case study demonstrating the effectiveness of joint multi-domain optimization. Finally, we discuss key challenges and future research directions to connect theoretical advancements and practical ISAC deployments. This work paves the way for intelligent and scalable ISAC architectures, providing critical insights for their seamless integration into next-generation wireless networks.

## 📚 Citation

If you use this simulation code package in any way, please cite our paper:

### BibTeX Format
```bibtex
@article{liu2025multidomain,
  title={Multi-domain optimization framework for {ISAC}: From electromagnetic shaping to network cooperation},
  author={Liu, Rang and Li, Ming and Zafari, Mehdi and Ottersten, Bj{\"o}rn and Swindlehurst, A. Lee},
  journal={IEEE Wireless Communications},
  year={2025},
  note={Under revision},
  url={https://arxiv.org/abs/2506.16011v1}
}
```

### IEEE Format
R. Liu, M. Li, M. Zafari, B. Ottersten, and A. L. Swindlehurst, "Multi-domain optimization framework for ISAC: From electromagnetic shaping to network cooperation," *IEEE Wireless Communications*, under revision, 2025. [Online]. Available: https://arxiv.org/abs/2506.16011v1

### Plain Text Format
Rang Liu, Ming Li, Mehdi Zafari, Björn Ottersten, and A. Lee Swindlehurst. "Multi-domain optimization framework for ISAC: From electromagnetic shaping to network cooperation." IEEE Wireless Communications (under revision), 2025.

## 🔧 Requirements

### Software Platform
- **MATLAB R2023b** or later (compatibility issues may arise with earlier versions)
- **[CVX](http://cvxr.com/cvx/)** - Convex optimization toolbox for solving SDP problems
- **[Manopt](https://www.manopt.org/)** - Optimization on manifolds toolbox for non-convex problems

### Installation
1. Download and install MATLAB R2023b or later
2. Install CVX toolbox
3. Install Manopt toolbox

## 📁 Code Structure

### Main Scripts
- **`main_power_SINR.m`** - Generates Figure 5(a): Radar SINR performance vs. transmit power
  - Sweeps transmit power from 30 to 50 dBm
  - Compares five optimization schemes
  - Monte Carlo simulation with 100 trials
  
- **`main_SR_SINR.m`** - Generates Figure 5(b): Radar SINR vs. communication sum-rate trade-off
  - Analyzes sensing-communication trade-offs
  - Fixed power at 50 dBm
  - Varies sum-rate requirement from 0 to 100 bps/Hz

### Optimization Functions

#### Core Optimization Algorithms
- **`opt_prop.m`** - Proposed cross-domain joint optimization
  - Jointly optimizes: AP selection, subcarrier allocation, beamforming, and polarization
  - Integrates EM, baseband, and network domains
  
- **`opt_radar.m`** - Radar-only benchmark
  - Maximizes radar SINR without communication constraints
  - Upper bound for sensing performance

#### Domain-Specific Optimization
- **`opt_BP.m`** - Baseband processing domain only
  - Fixed AP selection and polarization
  - Optimizes subcarrier allocation and beamforming
  
- **`opt_BP_EM.m`** - Joint baseband and EM-domain
  - Fixed AP selection
  - Optimizes subcarrier allocation, beamforming, and polarization
  
- **`opt_BP_NC.m`** - Joint baseband and network-domain
  - Fixed polarization
  - Optimizes AP selection, subcarrier allocation, and beamforming

### Channel Generation Functions
- **`generate_ap_to_ap_channels.m`** - Inter-AP interference channels (G_ijn)
- **`generate_ap_to_ue_channels.m`** - AP-to-UE communication channels (H_jnk)
- **`generate_target_scattering.m`** - Radar target scattering matrix (Φ_ij)

### Resource Allocation Functions
- **`allocate_ap.m`** - AP role assignment (transmit/receive selection)
- **`allocate_subcarrier.m`** - Frequency resource allocation between sensing and communication

### Compute Auxiliary Variables
- **`compute_c.m`**  
- **`compute_eta.m`**  
- **`compute_He.m`**  
- **`compute_t.m`**

### Initialization Functions
- **`initialize_MRT_grouped.m`** - Maximum ratio transmission initialization
- **`initialize_Wr_grouped.m`** - Radar beamforming initialization

### Update Functions  
- **`update_pk.m`** - User polarization
- **`update_rPj.m`** - Receive AP polarization
- **`update_tPj.m`** - Transmit AP polarization
- **`update_u.m`** - Receive filter 
- **`update_W.m`** - Transmit beamforming  

### Data and Visualization
- **`data_power.mat`** - Saved simulation results for power analysis
- **`data_sr.mat`** - Saved simulation results for sum-rate analysis
- **`power.fig`** / **`power.eps`** - Figure 5(a) outputs
- **`SR.fig`** / **`SR.eps`** - Figure 5(b) outputs

## 🎯 Key Features

### Multi-Domain Optimization
- **Electromagnetic Domain**: Polarization reconfigurability 
- **Baseband Domain**: Joint subcarrier allocation and beamforming design
- **Network Domain**: Collaborative AP selection 

### Performance Benchmarks
1. **Proposed**: Full cross-domain optimization
2. **Radar-only**: Sensing-optimal without communication
3. **BP-NC**: Baseband + Network cooperation
4. **BP-EM**: Baseband + EM shaping
5. **BP-only**: Baseband processing only


## 📈 Expected Results

- **Figure 5(a)**: Demonstrates up to 7 dB improvement in radar SINR through multi-domain optimization compared to baseband-only processing
- **Figure 5(b)**: Shows flexible sensing-communication trade-off with graceful performance degradation

## 👥 Authors

- **Rang Liu** (rangl2@uci.edu) - [Website](https://rangliu0706.github.io/)
  - Primary contributor and corresponding author
- **Ming Li** - Dalian University of Technology
- **Mehdi Zafari** - UC Irvine
- **Björn Ottersten** - University of Luxembourg
- **A. Lee Swindlehurst** - UC Irvine

## 📄 License

This code is licensed for personal, non-commercial use only, specifically for academic purposes. 

Copyright © 2025 LS Wireless Research Group  
Led by Prof. A. Lee Swindlehurst  
Department of Electrical Engineering and Computer Science  
University of California Irvine, CA 92697, USA

## 🔗 Related Resources

- **Paper**: [arXiv:2506.16011v1](https://arxiv.org/abs/2506.16011v1)
- **Research Group**: [LS-Wireless GitHub](https://github.com/LS-Wireless)
- **Author's Publications**: [https://rangliu0706.github.io/publications](https://rangliu0706.github.io/publications)

## 📞 Contact

For questions, suggestions, or collaboration inquiries:
- **Email**: rangl2@uci.edu

## 🙏 Acknowledgments

The authors would like to thank the reviewers for their valuable comments and suggestions. This work was supported in part by xxx.

---
*Last updated: September 28, 2025*
 