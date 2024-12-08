# Separator Sizing Calculator

A comprehensive GUI-based calculator for sizing separators in oil and gas operations. This open-source tool assists engineers and professionals in determining the appropriate separator sizes for various phase combinations, ensuring efficient and reliable separation processes.

## Features

- **Multi-Tab Interface**: 
  - **2-phase Vertical Sizing**
  - **2-phase Horizontal Sizing**
  - **3-phase Vertical Sizing**
  - **3-phase Horizontal Separator**
- **User-Friendly Inputs**: Easily enter and modify parameters through intuitive input fields.
- **Dynamic Calculations**: Perform complex sizing calculations based on industry-standard formulas.
- **Results Display**:
  - **Tables**: View calculated values for different diameters and configurations.
  - **Selected Separator Sizes**: Automatically highlights the most suitable separator sizes based on predefined criteria.
- **Error Handling**: Robust input validation and error messages to guide users.

## Installation

### Prerequisites

- **Python 3.x**: Ensure you have Python installed on your system. You can download it from [python.org](https://www.python.org/downloads/).

### Dependencies

The application relies on the following Python libraries:

- `tkinter`: Standard GUI library for Python (usually comes pre-installed).
- `ttk`: Themed widget set for `tkinter` (included with `tkinter`).
- `math`: Mathematical functions (standard library).

### Installation Steps

1. **Clone the Repository**

   ```bash
   git clone https://github.com/ViaAnthroposBenevolentia/GasOilSeparatorSizingCalculator.git
   cd GasOilSeparatorSizingCalculator
   ```

2. **(Optional) Create a Virtual Environment**

   It's good practice to use a virtual environment to manage dependencies.

   ```bash
   python3 -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   ```

3. **Install Dependencies**

   Since `tkinter` and `math` are part of the standard library, no additional installations are required. However, ensure that `tkinter` is installed on your system. If not, you can install it using:

   - **Ubuntu/Debian**

     ```bash
     sudo apt-get install python3-tk
     ```

   - **Fedora**

     ```bash
     sudo dnf install python3-tkinter
     ```

   - **macOS**

     `tkinter` comes pre-installed with Python on macOS.

## Usage

1. **Run the Application**

   Navigate to the project directory and execute the Python script:

   ```bash
   python gui_calculator.py
   ```

2. **Navigate Through Tabs**

   - **2-phase Vertical Sizing**: Calculate separator sizes for two-phase vertical systems.
   - **2-phase Horizontal Sizing**: Calculate separator sizes for two-phase horizontal systems.
   - **3-phase Vertical Sizing**: Calculate separator sizes for three-phase vertical systems.
   - **3-phase Horizontal Separator**: Calculate separator sizes for three-phase horizontal systems.

3. **Input Parameters**

   - Enter the required parameters in the input fields. Default values are provided based on standard examples.
   - Ensure all inputs are numerical to avoid calculation errors.

4. **Perform Calculations**

   - Click the **"Calculate"** button to execute the sizing calculations.
   - View the results in the displayed tables and the selected separator sizes based on the criteria.

## Contribution

Contributions are welcome! Whether it's reporting issues, suggesting features, or submitting pull requests, your input helps improve the project.

1. **Fork the Repository**

2. **Create a New Branch**

   ```bash
   git checkout -b feature/YourFeature
   ```

3. **Commit Your Changes**

   ```bash
   git commit -m "Add your message here"
   ```

4. **Push to the Branch**

   ```bash
   git push origin feature/YourFeature
   ```

5. **Create a Pull Request**

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgements

- Developed using Python's `tkinter` library.
- Inspired by industry-standard separator sizing methodologies.

## Contact

For any questions or suggestions, feel free to reach out:

- **Email**: alikhanov.2004@gmail.com