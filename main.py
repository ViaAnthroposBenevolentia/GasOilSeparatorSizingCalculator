import tkinter as tk
from tkinter import ttk, messagebox
import math

class SizingTab(ttk.Frame):
    def __init__(self, parent, labels, default_values, table_columns, table_headers):
        super().__init__(parent)
        self.labels = labels
        self.default_values = default_values
        self.table_columns = table_columns
        self.table_headers = table_headers
        self.entries = {}
        self.create_widgets()

    def create_widgets(self):
        # Input Frame
        input_frame = ttk.LabelFrame(self, text="Inputs")
        input_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

        # Create input fields
        for i, label in enumerate(self.labels):
            ttk.Label(input_frame, text=label).grid(row=i, column=0, sticky="w", padx=5, pady=2)
            entry = ttk.Entry(input_frame, validate="key", validatecommand=(self.register(self.validate_numeric), '%P'))
            entry.grid(row=i, column=1, padx=5, pady=2)
            entry.insert(0, self.default_values.get(label, ""))
            self.entries[label] = entry

        # Calculate Button
        calculate_button = ttk.Button(self, text="Calculate", command=self.calculate)
        calculate_button.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

        # Results Frame
        results_frame = ttk.LabelFrame(self, text="Results")
        results_frame.grid(row=2, column=0, padx=10, pady=10, sticky="nsew")

        # Configure grid to expand
        self.grid_rowconfigure(2, weight=1)
        self.grid_columnconfigure(0, weight=1)
        results_frame.grid_rowconfigure(0, weight=1)
        results_frame.grid_columnconfigure(0, weight=1)

        # Treeview for Table
        self.tree = ttk.Treeview(results_frame, columns=self.table_columns, show="headings")
        for col, header in zip(self.table_columns, self.table_headers):
            self.tree.heading(col, text=header)
            self.tree.column(col, anchor="center", width=100)
        self.tree.grid(row=0, column=0, sticky="nsew")

        # Scrollbar for Treeview
        scrollbar = ttk.Scrollbar(results_frame, orient=tk.VERTICAL, command=self.tree.yview)
        self.tree.configure(yscroll=scrollbar.set)
        scrollbar.grid(row=0, column=1, sticky='ns')

        # Result Text
        self.result_text = tk.Text(results_frame, height=6, state='disabled', wrap='word')
        self.result_text.grid(row=1, column=0, columnspan=2, padx=5, pady=5, sticky="ew")

    def validate_numeric(self, P):
        """Validate that the input is a valid float or empty."""
        if P == "":
            return True
        try:
            float(P)
            return True
        except ValueError:
            return False

    def calculate(self):
        """To be implemented by subclasses."""
        raise NotImplementedError("Calculate method must be implemented by the subclass.")

    def clear_treeview(self):
        """Clear all items from the treeview."""
        for item in self.tree.get_children():
            self.tree.delete(item)

    def insert_into_treeview(self, data):
        """Insert data into the treeview."""
        for row in data:
            values = tuple(row[col] for col in self.table_columns)
            self.tree.insert('', 'end', values=values)

    def display_result(self, text):
        """Display result text in the result_text widget."""
        self.result_text.config(state='normal')
        self.result_text.delete(1.0, tk.END)
        self.result_text.insert(tk.END, text)
        self.result_text.config(state='disabled')

    def iterative_cd_calculation(self, rho_l, rho_g, d_m, mu):
        """Iteratively calculate C_D."""
        C_D = 0.34
        tolerance = 1e-3
        max_iterations = 100
        for _ in range(max_iterations):
            term = ((rho_l - rho_g) / rho_g) * (d_m / C_D)
            V_t = 0.0119 * math.sqrt(term)
            Re = 0.0049 * ((rho_g * d_m * V_t) / mu)
            if Re == 0:
                raise ZeroDivisionError("Reynolds number calculation resulted in zero.")
            C_D_new = 24 / Re + 3 / math.sqrt(Re) + 0.34
            if abs(C_D_new - C_D) < tolerance:
                return C_D_new
            C_D = C_D_new
        messagebox.showwarning("Calculation Warning", "C_D did not converge within the maximum number of iterations.")
        return None

class TwoPhaseVerticalSizingTab(SizingTab):
    def __init__(self, parent):
        labels = [
            "Gas Flow (MMSCFD):",
            "Specific Gravity (Gas):",
            "Oil Flow (BOPD):",
            "API Gravity:",
            "Pressure (psia):",
            "Temperature (°F):",
            "Particle Size (microns):",
            "Retention Time (minutes):",
            "TZ:",
            "Z Factor:",
            "Viscosity (cp):"
        ]

        default_values = {
            "Gas Flow (MMSCFD):": "10",
            "Specific Gravity (Gas):": "0.6",
            "Oil Flow (BOPD):": "2000",
            "API Gravity:": "40",
            "Pressure (psia):": "1000",
            "Temperature (°F):": "60",
            "Particle Size (microns):": "140",
            "Retention Time (minutes):": "3",
            "TZ:": "520",
            "Z Factor:": "0.84",
            "Viscosity (cp):": "0.013"
        }

        table_columns = ("t_r (min)", "d (in.)", "h (in.)", "L_ss (ft.)", "S_R")
        table_headers = ("Retention Time (min)", "Diameter (in.)", "Height (in.)", "Seam-to-Seam Length (ft.)", "Slenderness Ratio (S_R)")

        super().__init__(parent, labels, default_values, table_columns, table_headers)

    def calculate(self):
        try:
            # Retrieve and parse inputs
            inputs = {label: float(entry.get()) for label, entry in self.entries.items()}

            gas_flow_mmscfd = inputs["Gas Flow (MMSCFD):"]
            gas_sp = inputs["Specific Gravity (Gas):"]
            oil_flow_bopd = inputs["Oil Flow (BOPD):"]
            api = inputs["API Gravity:"]
            P = inputs["Pressure (psia):"]
            T = inputs["Temperature (°F):"]
            d_m = inputs["Particle Size (microns):"]
            retention_time = inputs["Retention Time (minutes):"]
            TZ = inputs["TZ:"]
            Z = inputs["Z Factor:"]
            mu = inputs["Viscosity (cp):"]

            # Calculate liquid and gas densities
            rho_l = 62.4 * (141.5 / (131.5 + api))  # lb/ft³
            rho_g = 2.70 * (gas_sp * P) / (TZ * Z)  # lb/ft³

            # Iterative calculation for C_D
            C_D = self.iterative_cd_calculation(rho_l, rho_g, d_m, mu)
            if C_D is None:
                return

            # Gas capacity constraint
            d_squared = 5.040 * ((TZ * gas_flow_mmscfd) / P) * (((rho_l - rho_g) * C_D) / d_m) ** 0.5
            d_gas = math.sqrt(d_squared)  # inches

            # Liquid capacity constraint
            t_r_list = [3, 2, 1]  # minutes
            d_list = [24, 30, 36, 42, 48]  # inches

            # Prepare table data
            table = []
            for t_r in t_r_list:
                for d in d_list:
                    h = (t_r * oil_flow_bopd) / (0.12 * d ** 2)  # inches
                    L_ss = (h + 76) / 12  # feet
                    S_R = (12 * L_ss) / d
                    table.append({
                        't_r (min)': t_r,
                        'd (in.)': d,
                        'h (in.)': round(h, 1),
                        'L_ss (ft.)': round(L_ss, 1),
                        'S_R': round(S_R, 1)
                    })

            # Update treeview
            self.clear_treeview()
            self.insert_into_treeview(table)

            # Gas capacity diameter
            gas_capacity_text = f"Gas capacity diameter: {round(d_gas, 1)} in."

            # Select suitable separator size
            selected = next((row for row in table if row['d (in.)'] >= d_gas and 3 <= row['S_R'] <= 4), None)

            if selected:
                selected_text = (
                    f"\nSelected Separator Size:\n"
                    f"Retention Time: {selected['t_r (min)']} minutes\n"
                    f"Diameter: {selected['d (in.)']} in.\n"
                    f"Seam-to-Seam Length: {selected['L_ss (ft.)']} ft.\n"
                    f"Slenderness Ratio: {selected['S_R']}"
                )
            else:
                selected_text = "\nNo suitable separator size found."

            # Display results
            self.display_result(gas_capacity_text + selected_text)

        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numerical values.")
        except ZeroDivisionError as e:
            messagebox.showerror("Calculation Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

class TwoPhaseHorizontalSizingTab(SizingTab):
    def __init__(self, parent):
        labels = [
            "Gas Flow (MMSCFD):",
            "Specific Gravity (Gas):",
            "Oil Flow (BOPD):",
            "API Gravity:",
            "Operating Pressure (psia):",
            "Operating Temperature (°F):",
            "Particle Size (microns):",
            "Retention Time (minutes):",
            "T:",
            "Z Factor:",
            "Viscosity (cp):"
        ]

        default_values = {
            "Gas Flow (MMSCFD):": "10",
            "Specific Gravity (Gas):": "0.6",
            "Oil Flow (BOPD):": "2000",
            "API Gravity:": "40",
            "Operating Pressure (psia):": "1000",
            "Operating Temperature (°F):": "60",
            "Particle Size (microns):": "140",
            "Retention Time (minutes):": "3",
            "T:": "520",
            "Z Factor:": "0.84",
            "Viscosity (cp):": "0.013"
        }

        table_columns = ("d (in.)", "Gas L_eff (ft)", "Liquid L_eff (ft)", "L_ss (ft.)", "S_R")
        table_headers = ("Diameter (in.)", "Gas L_eff (ft)", "Liquid L_eff (ft)", "Seam-to-Seam Length (ft.)", "Slenderness Ratio (S_R)")

        super().__init__(parent, labels, default_values, table_columns, table_headers)

    def calculate(self):
        try:
            # Retrieve and parse inputs
            inputs = {label: float(entry.get()) for label, entry in self.entries.items()}

            gas_flow_mmscfd = inputs["Gas Flow (MMSCFD):"]
            gas_sp = inputs["Specific Gravity (Gas):"]
            oil_flow_bopd = inputs["Oil Flow (BOPD):"]
            api = inputs["API Gravity:"]
            O_P = inputs["Operating Pressure (psia):"]
            O_T = inputs["Operating Temperature (°F):"]
            d_m = inputs["Particle Size (microns):"]
            retention_time = inputs["Retention Time (minutes):"]
            T = inputs["T:"]
            Z = inputs["Z Factor:"]
            mu = inputs["Viscosity (cp):"]

            # Calculate liquid and gas densities
            rho_l = 62.4 * (141.5 / (131.5 + api))  # lb/ft³
            rho_g = 2.70 * (gas_sp * O_P) / (T * Z)  # lb/ft³

            # Iterative calculation for C_D
            C_D = self.iterative_cd_calculation(rho_l, rho_g, d_m, mu)
            if C_D is None:
                return

            # Gas capacity constraint
            dL_eff = 420 * ((T * Z * gas_flow_mmscfd) / O_P) * math.sqrt((rho_g / (rho_l - rho_g)) * (C_D / d_m))

            # Prepare table data
            d_list = [16, 20, 24, 30, 36, 42, 48]  # in inches
            Q_l = oil_flow_bopd
            t_r = retention_time
            table = []
            for d in d_list:
                L_eff_gas = dL_eff / d  # ft
                L_eff_liquid = (t_r * Q_l) / (0.7 * d ** 2)  # ft
                L_ss = L_eff_liquid + (d / 12)  # ft
                S_R = (12 * L_ss) / d
                table.append({
                    'd (in.)': d,
                    'Gas L_eff (ft)': round(L_eff_gas, 1),
                    'Liquid L_eff (ft)': round(L_eff_liquid, 1),
                    'L_ss (ft.)': round(L_ss, 1),
                    'S_R': round(S_R, 1)
                })

            # Update treeview
            self.clear_treeview()
            self.insert_into_treeview(table)

            # Select suitable separator size
            selected = next((row for row in table if 3 <= row['S_R'] <= 4 and row['d (in.)'] >= dL_eff ** 0.5), None)

            if selected:
                selected_text = (
                    f"\nSelected Separator Size:\n"
                    f"Diameter: {selected['d (in.)']} in.\n"
                    f"Seam-to-Seam Length: {selected['L_ss (ft.)']} ft.\n"
                    f"Gas L_eff: {selected['Gas L_eff (ft)']} ft.\n"
                    f"Liquid L_eff: {selected['Liquid L_eff (ft)']} ft.\n"
                    f"Slenderness Ratio: {selected['S_R']}"
                )
            else:
                selected_text = "\nNo suitable separator size found."

            # Display results
            self.display_result(selected_text)

        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numerical values.")
        except ZeroDivisionError as e:
            messagebox.showerror("Calculation Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

class ThreePhaseVerticalSizingTab(SizingTab):
    def __init__(self, parent):
        labels = [
            "Oil Flow (BOPD):",
            "Water Flow (BWPD):",
            "Gas Flow (MMscfd):",
            "Operating Pressure (psia):",
            "Operating Temperature (°F):",
            "API Gravity:",
            "Specific Gravity of Water:",
            "Specific Gravity of Gas:",
            "Retention Time Oil (min):",
            "Retention Time Water (min):",
            "Viscosity Oil (cp):",
            "Viscosity Water (cp):",
            "Gas Density (lb/ft³):",
            "Liquid Density (lb/ft³):",
            "C_D:",
            "Droplet Removal Liquid (microns):",
            "Droplet Removal Water (microns):",
            "Droplet Removal Oil (microns):"
        ]

        default_values = {
            "Oil Flow (BOPD):": "5000",
            "Water Flow (BWPD):": "3000",
            "Gas Flow (MMscfd):": "5",
            "Operating Pressure (psia):": "100",
            "Operating Temperature (°F):": "90",
            "API Gravity:": "30",
            "Specific Gravity of Water:": "1.07",
            "Specific Gravity of Gas:": "0.99",
            "Retention Time Oil (min):": "10",
            "Retention Time Water (min):": "10",
            "Viscosity Oil (cp):": "10",
            "Viscosity Water (cp):": "1",
            "Gas Density (lb/ft³):": "0.3",
            "Liquid Density (lb/ft³):": "54.7",
            "C_D:": "2.01",
            "Droplet Removal Liquid (microns):": "100",
            "Droplet Removal Water (microns):": "500",
            "Droplet Removal Oil (microns):": "200"
        }

        table_columns = ("d_o (in.)", "h_o + h_w (in.)", "L_ss (ft)", "SR (12L_ss / d_o)")
        table_headers = ("Diameter (in.)", "Height (in.)", "Seam-to-Seam Length (ft)", "Slenderness Ratio (S_R)")

        super().__init__(parent, labels, default_values, table_columns, table_headers)

    def calculate(self):
        try:
            # Retrieve and parse inputs
            inputs = {label: float(entry.get()) for label, entry in self.entries.items()}

            Q_o = inputs["Oil Flow (BOPD):"]
            Q_w = inputs["Water Flow (BWPD):"]
            Q_g = inputs["Gas Flow (MMscfd):"]
            P_o = inputs["Operating Pressure (psia):"]
            T_o = inputs["Operating Temperature (°F):"]
            API = inputs["API Gravity:"]
            SG_w = inputs["Specific Gravity of Water:"]
            S_g = inputs["Specific Gravity of Gas:"]
            t_r_o = inputs["Retention Time Oil (min):"]
            t_r_w = inputs["Retention Time Water (min):"]
            mu_o = inputs["Viscosity Oil (cp):"]
            mu_w = inputs["Viscosity Water (cp):"]
            rho_g = inputs["Gas Density (lb/ft³):"]
            rho_l = inputs["Liquid Density (lb/ft³):"]
            C_D = inputs["C_D:"]
            droplet_removal_liquid = inputs["Droplet Removal Liquid (microns):"]
            droplet_removal_water = inputs["Droplet Removal Water (microns):"]
            droplet_removal_oil = inputs["Droplet Removal Oil (microns):"]

            # Step 1: Calculate difference in specific gravities
            SG_o = 141.5 / (API + 131.5)  # Specific Gravity of oil
            Delta_SG = SG_w - SG_o
            if Delta_SG == 0:
                raise ZeroDivisionError("Delta_SG is zero, cannot proceed with calculations.")

            # Step 2: Calculate minimum diameter for liquid droplet settling through gas phase
            TZ = 550  # Assumed constant
            d_m_liquid = droplet_removal_liquid  # microns
            d2_liquid_settle = 5040 * ((TZ * S_g * Q_g) / P_o) * math.sqrt((rho_g / (rho_l - rho_g)) * (C_D / d_m_liquid))
            if d2_liquid_settle < 0:
                raise ValueError("Calculated d2_liquid_settle is negative. Check input values.")
            d_liquid_settle = math.sqrt(d2_liquid_settle)  # inches

            # Step 3: Calculate minimum diameter for water droplets to settle through oil phase
            d_m_water = droplet_removal_water  # microns
            d2_water_settle = 6690 * ((Q_o * mu_o) / (Delta_SG * (d_m_water ** 2)))
            if d2_water_settle < 0:
                raise ValueError("Calculated d2_water_settle is negative. Check input values.")
            d_water_settle = math.sqrt(d2_water_settle)  # inches

            # Step 4: Calculate minimum diameter for oil droplets to rise through water phase
            d_m_oil = droplet_removal_oil  # microns
            d2_oil_rise = 6690 * ((Q_w * mu_w) / (Delta_SG * (d_m_oil ** 2)))
            if d2_oil_rise < 0:
                raise ValueError("Calculated d2_oil_rise is negative. Check input values.")
            d_oil_rise = math.sqrt(d2_oil_rise)  # inches

            # Step 5: Select the largest diameter as the minimum inside diameter required
            d_min = max(d_liquid_settle, d_water_settle, d_oil_rise)

            # Step 6: Calculate h_o + h_w
            h_o_w_constant = (t_r_o * Q_o + t_r_w * Q_w) / 0.12  # inches

            # Define a list of diameters to evaluate (in inches)
            diameters = [84, 90, 96, 102]  # in inches

            # Prepare table data
            table = []
            for d in diameters:
                h_o_w = h_o_w_constant / (d ** 2)  # inches
                # Compute seam-to-seam length (L_ss)
                if d <= 36:
                    L_ss = (h_o_w + 76) / 12  # feet
                else:
                    L_ss = (h_o_w + d + 40) / 12  # feet
                # Compute slenderness ratio
                S_R = (12 * L_ss) / d
                table.append({
                    'd_o (in.)': d,
                    'h_o + h_w (in.)': round(h_o_w, 1),
                    'L_ss (ft)': round(L_ss, 1),
                    'SR (12L_ss / d_o)': round(S_R, 1)
                })

            # Update treeview
            self.clear_treeview()
            self.insert_into_treeview(table)

            # Step 7: Make final selection
            selected = next((row for row in table if row['d_o (in.)'] >= d_min and 1.5 <= row['SR (12L_ss / d_o)'] <= 3), None)

            if selected:
                selected_text = (
                    f"Selected Separator Size:\n"
                    f"Diameter: {selected['d_o (in.)']} in.\n"
                    f"Seam-to-seam Length: {selected['L_ss (ft)']} ft.\n"
                    f"Retention Time Oil: {t_r_o} minutes\n"
                    f"Slenderness Ratio: {selected['SR (12L_ss / d_o)']}"
                )
            else:
                selected_text = "No suitable separator size found."

            # Display results
            self.display_result(selected_text)

        except ValueError as ve:
            messagebox.showerror("Calculation Error", f"Value Error: {ve}")
        except ZeroDivisionError as zde:
            messagebox.showerror("Calculation Error", f"Zero Division Error: {zde}")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

class ThreePhaseHorizontalSizingTab(SizingTab):
    def __init__(self, parent):
        labels = [
            "Oil Flow (BOPD):",
            "Water Flow (BWPD):",
            "Gas Flow (MMscfd):",
            "Operating Pressure (psia):",
            "Operating Temperature (°F):",
            "API Gravity:",
            "Specific Gravity of Water:",
            "Specific Gravity of Gas:",
            "Retention Time Oil (min):",
            "Retention Time Water (min):",
            "Viscosity Oil (cp):",
            "Viscosity Water (cp):",
            "Gas Density (lb/ft³):",
            "Liquid Density (lb/ft³):",
            "C_D:",
            "Droplet Removal Liquid (microns):",
            "Droplet Removal Water (microns):",
            "Droplet Removal Oil (microns):"
        ]

        default_values = {
            "Oil Flow (BOPD):": "5000",
            "Water Flow (BWPD):": "3000",
            "Gas Flow (MMscfd):": "5",
            "Operating Pressure (psia):": "100",
            "Operating Temperature (°F):": "90",
            "API Gravity:": "30",
            "Specific Gravity of Water:": "1.07",
            "Specific Gravity of Gas:": "0.99",
            "Retention Time Oil (min):": "10",
            "Retention Time Water (min):": "10",
            "Viscosity Oil (cp):": "10",
            "Viscosity Water (cp):": "1",
            "Gas Density (lb/ft³):": "0.3",
            "Liquid Density (lb/ft³):": "54.7",
            "C_D:": "2.01",
            "Droplet Removal Liquid (microns):": "100",
            "Droplet Removal Water (microns):": "500",
            "Droplet Removal Oil (microns):": "200"
        }

        table_columns = ("d_o (in.)", "h_o + h_w (in.)", "L_ss (ft)", "SR (12L_ss / d_o)")
        table_headers = ("Diameter (in.)", "Height (in.)", "Seam-to-Seam Length (ft)", "Slenderness Ratio (S_R)")

        super().__init__(parent, labels, default_values, table_columns, table_headers)

    def interpolate_beta(self, A_w_over_A, A_list, beta_list):
        """Interpolate beta based on A_w_over_A using linear interpolation."""
        if A_w_over_A <= A_list[0]:
            # Extrapolate below the first point
            slope = (beta_list[1] - beta_list[0]) / (A_list[1] - A_list[0])
            beta = beta_list[0] + slope * (A_w_over_A - A_list[0])
        elif A_w_over_A >= A_list[-1]:
            # Extrapolate above the last point
            slope = (beta_list[-1] - beta_list[-2]) / (A_list[-1] - A_list[-2])
            beta = beta_list[-1] + slope * (A_w_over_A - A_list[-1])
        else:
            # Interpolate between points
            for i in range(len(A_list) - 1):
                if A_list[i] <= A_w_over_A <= A_list[i + 1]:
                    slope = (beta_list[i + 1] - beta_list[i]) / (A_list[i + 1] - A_list[i])
                    beta = beta_list[i] + slope * (A_w_over_A - A_list[i])
                    break
        return beta

    def calculate(self):
        try:
            # Retrieve and parse inputs
            inputs = {label: float(entry.get()) for label, entry in self.entries.items()}

            Q_o = inputs["Oil Flow (BOPD):"]
            Q_w = inputs["Water Flow (BWPD):"]
            Q_g = inputs["Gas Flow (MMscfd):"]
            P_s = inputs["Operating Pressure (psia):"]
            T = inputs["Operating Temperature (°F):"]
            API = inputs["API Gravity:"]
            S_w = inputs["Specific Gravity of Water:"]
            S_g = inputs["Specific Gravity of Gas:"]
            t_r_o = inputs["Retention Time Oil (min):"]
            t_r_w = inputs["Retention Time Water (min):"]
            mu_o = inputs["Viscosity Oil (cp):"]
            mu_w = inputs["Viscosity Water (cp):"]
            rho_g = inputs["Gas Density (lb/ft³):"]
            rho_l = inputs["Liquid Density (lb/ft³):"]
            C_D = inputs["C_D:"]
            droplet_removal_liquid = inputs["Droplet Removal Liquid (microns):"]
            droplet_removal_water = inputs["Droplet Removal Water (microns):"]
            droplet_removal_oil = inputs["Droplet Removal Oil (microns):"]

            # Step 1: Calculate specific gravity of oil and ΔSG
            SG_o = 141.5 / (API + 131.5)  # Specific Gravity of oil
            Delta_SG = S_w - SG_o
            if Delta_SG == 0:
                raise ZeroDivisionError("Delta_SG is zero, cannot proceed with calculations.")

            # Step 2: Calculate maximum oil pad thickness (h_o_max)
            h_o_max = 0.00128 * t_r_o * Delta_SG * (droplet_removal_water ** 2) / mu_o  # in inches

            # Step 3: Calculate A_w / A
            A_w_over_A = 0.5 * Q_w / (Q_o + Q_w)

            # Step 4: Determine β from Figure 5-20 (interpolation)
            A_list = [0.1, 0.2, 0.3]
            beta_list = [0.3, 0.25, 0.2]
            beta = self.interpolate_beta(A_w_over_A, A_list, beta_list)

            # Step 5: Calculate d_max
            d_max = h_o_max / beta  # in inches

            # Step 6: Gas capacity constraint (Skipped as per example, gas does not govern)

            # Step 7: Liquid retention time constraint
            constant = 1.42 * (Q_o * t_r_o + Q_w * t_r_w)

            # Step 8: Estimate seam-to-seam length
            d_values = [60, 72, 84, 96, 108]  # in inches
            table = []
            for d in d_values:
                L_eff = constant / (d ** 2)  # in feet
                L_ss = (4 * L_eff) / 3  # in feet
                SR = (12 * L_ss) / d
                table.append({
                    'd_o (in.)': d,
                    'h_o + h_w (in.)': round(L_eff * 0.3, 1),  # Assuming h_o + h_w relates to L_eff
                    'L_ss (ft)': round(L_ss, 1),
                    'SR (12L_ss / d_o)': round(SR, 1)
                })

            # Update treeview
            self.clear_treeview()
            self.insert_into_treeview(table)

            # Step 9: Make final selection
            selected_sizes = [
                row for row in table
                if 3 <= row['SR (12L_ss / d_o)'] <= 5 and row['d_o (in.)'] <= d_max
            ]

            # Step 10: Choose reasonable size
            if selected_sizes:
                selected_text = "Selected Separator Sizes:\n"
                for size in selected_sizes:
                    selected_text += (
                        f"Diameter: {size['d_o (in.)']} in.\n"
                        f"Seam-to-seam Length: {size['L_ss (ft)']} ft.\n"
                        f"Retention Time Oil: {t_r_o} minutes\n"
                        f"Slenderness Ratio: {size['SR (12L_ss / d_o)']}\n\n"
                    )
            else:
                selected_text = "No suitable separator size found."

            # Display results
            self.display_result(selected_text)

        except ValueError as ve:
            messagebox.showerror("Calculation Error", f"Value Error: {ve}")
        except ZeroDivisionError as zde:
            messagebox.showerror("Calculation Error", f"Zero Division Error: {zde}")
        except Exception as e:
            messagebox.showerror("Error", f"An unexpected error occurred: {e}")

class CalculatorApp(ttk.Frame):
    def __init__(self, parent):
        super().__init__(parent)
        self.create_widgets()

    def create_widgets(self):
        notebook = ttk.Notebook(self)
        notebook.pack(expand=1, fill="both")

        # 2-phase Vertical Sizing Tab
        two_phase_vertical_tab = TwoPhaseVerticalSizingTab(notebook)
        notebook.add(two_phase_vertical_tab, text="2-phase Vertical Sizing")

        # 2-phase Horizontal Sizing Tab
        two_phase_horizontal_tab = TwoPhaseHorizontalSizingTab(notebook)
        notebook.add(two_phase_horizontal_tab, text="2-phase Horizontal Sizing")

        # 3-phase Vertical Sizing Tab
        three_phase_vertical_tab = ThreePhaseVerticalSizingTab(notebook)
        notebook.add(three_phase_vertical_tab, text="3-phase Vertical Sizing")

        # 3-phase Horizontal Sizing Tab
        three_phase_horizontal_tab = ThreePhaseHorizontalSizingTab(notebook)
        notebook.add(three_phase_horizontal_tab, text="3-phase Horizontal Sizing")

def main():
    root = tk.Tk()
    root.title("Separator Sizing Calculator")
    root.geometry("1200x800")
    app = CalculatorApp(root)
    app.pack(expand=True, fill="both")
    root.mainloop()

if __name__ == "__main__":
    main()