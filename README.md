# Apex Slicer

A 4-axis non-planar slicer, optimized for execution on Jupyter Notebook.

This project is designed for users who want a more streamlined execution and simplified workflow for 4-axis printing.

## 🚀 Key Features
- **Enhanced Compatibility:** Optimized for high-end GPUs (tested on RTX 3080 Ti) to handle complex meshes and "Torture Tests" without kernel crashes.
- **Progress Bars:** Added progress bars so you don't think nothing is happening (it's just extremely slow).
- **Standard Commands:** Works with `G01 Cxx Xxx Zxx Bxx Fxx` G-code structure.
- **Porting to GPU:** Some parts of Apex now works with GPU instead of CPU, ensuring stunning operation speeds (16 mins to 30 secs).
- **In Development:** A standalone UI EXE version is currently being cooked! 👨‍💻

## 📺 Media & Tutorials
Check out the project in action and learn how to use it:

[![Watch the video](https://github.com/aeloxhardware/apex/blob/main/thumbnail.webp?raw=true)](https://aeloxhardware.com)

## 🤝 Community & Support
Need help or want to share your results?
Please use the [Discussions tab](https://github.com/aeloxhardware/apex/discussions) to ask questions and interact with the community.

## 📈 Development Roadmap & Metrics
**Standalone EXE Interface Pre-Release:** ![Progress](https://geps.dev/progress/95)

**Source Code Adoption Goal:** ![Clones](https://img.shields.io/badge/Goal-12%20%2F%2050%20clones-blue?style=for-the-badge&logo=github&color=blue)

**Total Releases Downloads:**![Total Downloads](https://img.shields.io/github/downloads/aeloxhardware/apex/total?style=for-the-badge&logo=github&color=31c653)

## Run
You should have Python 3.12.x and CUDA 12.x. Go to the root folder, open PowerShell INTO THIS FOLDER and type 'python run main.py'. If dependencies are missing, install them.

### Bibtex Citation:
```bibtex
@software{S4_Slicer},
author = {Bird, Joshua},
license = {GPL-3.0},
title = {{S4 Slicer}},
url = {[https://github.com/jyjblrd/S4_Slicer](https://github.com/jyjblrd/S4_Slicer)}
}
