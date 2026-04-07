import streamlit as st
from CoolProp.CoolProp import PropsSI
from numpy import *
from scipy.optimize import root_scalar
from matplotlib.pyplot import *
import uuid
import pandas as pd

st.sidebar.write("Some Default Values")
Kelvin=st.sidebar.text_input("Kelvin conversion",value='273.15')
Kelvin = float(Kelvin)

Hf_Water = st.sidebar.text_input("Enthalpy (△H) of Water Formation at Standard Condition (in J/mol)",value='-285840')
Hf_Water = float(Hf_Water)    #J/mol

Gf_Water = st.sidebar.text_input("Gibbs Free Energy(△G) of Water Formation at Standard Condition (in J/mol)",value='-237220')
Gf_Water = float(Gf_Water)    #J/mol

St_T = st.sidebar.text_input("Temperature at Standard Condition (°C)",value='25')
St_T = float(St_T)     # in degree C

St_P = st.sidebar.text_input("Pressure at Standard Condition (bar)",value='1')
St_P = float(St_P)      # in bar

F = st.sidebar.text_input("Faraday's Constant (C/mol)",value='96485')
F = float(F)     # C/mol

n = st.sidebar.text_input("Number of electrons transferred",value='2')
n = float(n)         # number of electron transfer

M_H2O = st.sidebar.text_input("Molecular Weight of Water (g/mol)",value='18')
M_H2O = float(M_H2O)    # g/mol

R = st.sidebar.text_input("Ideal Gas Constant (J/(mol.K))",value='8.314')
R = float(R)             # Ideal gas constant in J/mol/K


# Operational Parameters (determined or set before experiment)

st.title("Single Cell Electrolysis Simulation")
st.divider()
st.write("This page simulates a single electrolysis cell, small enough to ignore the temperature inlet and outlet difference.")

st.markdown("**How to configure?**")
st.write("First, you need to enter the parameters below. After filling these, the results will be shown in two parts: single value results, and the plots. You may change a single parameter by checking **Run multi-factor experiment** and the new results will be shown in the same table or plot by adding values to that parameter. You have to configure your plot settings to get the desirable plot(s). Multiple y-axes can be selected and **multi-factor** could also be run (It is advised to not do both in the same plot). Some default values are given in the sidebar on the left. You may change them if required.")


st.header("Enter the operating parameters")


#st.write

with st.container(horizontal=True):
    st.write("Operating Temperature:")
    T_op = st.text_input("Temperature",placeholder = "in °C",label_visibility="collapsed")   # in degree C
T_op = float(T_op) if T_op else 0
T_op_K = T_op + Kelvin     # in Kelvin
st.caption("This temperature is the operating temperature and it will be assumed that there is negligible temperature difference in inlet and outlet of the cell.")
with st.container(horizontal=True):
    st.write("Operating Pressure:")
    P_op = st.text_input("Pressure",placeholder = "in bar",label_visibility="collapsed")
P_op = float(P_op) if P_op else 1
st.caption("The pressure is applied to both ends of the cell i.e. anode and cathode.")
with st.container(horizontal=True):
    st.write("Final Operating Voltage:")
    V_final = st.text_input("Voltage",placeholder = "in Volts",label_visibility="collapsed")
V_final = float(V_final) if V_final else 0   # Final operating voltage
st.caption("In the model, the voltage will be ramped up to this final operating voltage from reversible voltage.")
with st.container(horizontal=True):
    st.write("Limiting Current Density:")
    j_lim = st.text_input("j_lim",placeholder = "in Ampere/cm²",label_visibility="collapsed")
j_lim = float(j_lim) if j_lim else 0           # limiting current density in A/cm2
st.caption("This is the maximum current density that can be achived.")
with st.container(horizontal=True):
    st.write("Mass Activity of IrO₂:")
    Mass_activity_IrO2 = st.text_input("Mass activity",placeholder = "in Ampere/mg",label_visibility="collapsed")
st.caption("In a PEM electrolyser, the oxygen evolution reaction is the rate determing step, hence its mass activity is required to be known. You may write 0.0000015 as 1.5e-6")
Mass_activity_IrO2 = float(Mass_activity_IrO2) if Mass_activity_IrO2 else 0      # Ampere/mg

with st.container(horizontal=True):
    st.write("IrO₂ Loading:")
    Loading_IrO2 = st.text_input("Loading",placeholder = "in mg/cm² (Usually 1 to 3)",label_visibility="collapsed")
Loading_IrO2 =  float(Loading_IrO2) if Loading_IrO2 else 0      # mg/cm2 catalyst loading
st.caption("Catalyst loading is a common term in the field of catalysis. It tells how much catalyst has been added per sq. cms. of the cell area.")
#j0_ref = 1e-6     # exchange current density in A/cm2. Usually in range 1e-7 to 1e-9
with st.container(horizontal=True):
    st.write("IrO₂ Activation Energy at 25°C:")
    j0_ref_E_act = st.text_input("Activation Energy",placeholder = "in kJ/mol",label_visibility="collapsed")
st.caption("Activation energy is the minimum energy required to start the reaction. Usually in the range 60 to 80 kJ/mol")
j0_ref_E_act = float(j0_ref_E_act) if j0_ref_E_act else 0  # kJ/mol. Usually in the range 60 to 80 kJ/mol
#j0_ref_E_act = j0_ref_E_act*1000
with st.container(horizontal=True):
    st.write("Asymmetric Factor (ß)")
    beta = st.text_input("beta",placeholder = "",label_visibility="collapsed")
st.caption("ranges from 0 to 1, usually 0.5 is used to show symmetry")
beta = float(beta) if beta else 0           # asymmetric factor

with st.container(horizontal=True):
    st.write("Membrane Weight%:")
    wt_percent_membrane = st.text_input("mem wt percent",placeholder = "",label_visibility="collapsed")
st.caption("It is ratio water hold up of membrane to its dry weight. It ranges from 0.3 to 0.38 for fully hydrated membrane in 80-100°C")
wt_percent_membrane = float(wt_percent_membrane) if wt_percent_membrane else 0     # usually in the range 0.3 to 0.38 for a fully hydrated membrane near boiling conditions.

with st.container(horizontal=True):
    st.write("Membrane Equivalent Weight(EW):")
    EW = st.text_input("EW",placeholder = "g/mol",label_visibility="collapsed")
st.caption("Equivalent weight of nafion 117 membrane. Defined as membrane weight in gram per mol of -SO3H. Usually in the range 800-1500.")
EW = float(EW) if EW else 0            # Equivalent weight of nafion 117 membrane. Defined as membrane weight in gram per mol of -SO3H. Usually in the range 800-1500.

with st.container(horizontal=True):
    st.write("Membrane Thickness:")
    t_membrane = st.text_input("mem thickness",placeholder = "in microns (write '200' if thickness is 200 microns)",label_visibility="collapsed") 
t_membrane = float(t_membrane) if t_membrane else 0       # membrane thickness in microns (1e-6m)
st.caption("It is the hydrated thickness of membrane")

def Electrolysis_simulation(T_op,P_op,V_final,j_lim,Mass_activity_IrO2,Loading_IrO2,j0_ref_E_act,beta,wt_percent_membrane,EW,t_membrane):
    j0_ref_E_act = j0_ref_E_act*1000
    def get_enthalpy(fluid,T,P):
        return PropsSI('Hmolar','T',(T+Kelvin),'P',P*1e5,fluid)

    def get_gibbs(fluid,T,P):
        return PropsSI('Gmolar','T',(T+Kelvin),'P',P*1e5,fluid)

    del_H_H2 = get_enthalpy('Hydrogen',T_op,P_op) - get_enthalpy('Hydrogen',St_T,St_P)
    del_H_O2 = get_enthalpy('Oxygen',T_op,P_op) - get_enthalpy('Oxygen',St_T,St_P)
    del_H_H2O = Hf_Water + get_enthalpy('Water',T_op,P_op) - get_enthalpy('Water',St_T,St_P)

    Reaction_del_H = del_H_H2 + 0.5*del_H_O2 - del_H_H2O
    Thermoneutral_Voltage = Reaction_del_H/(n*F)

    del_G_H2 = get_gibbs('Hydrogen',T_op,P_op) - get_gibbs('Hydrogen',St_T,St_P)
    del_G_O2 = get_gibbs('Oxygen',T_op,P_op) - get_gibbs('Oxygen',St_T,St_P)
    del_G_H2O = Gf_Water + get_gibbs('Water',T_op,P_op) - get_gibbs('Water',St_T,St_P)

    Reaction_del_G = del_G_H2 + 0.5*del_G_O2 - del_G_H2O
    Reversible_Voltage = Reaction_del_G/(n*F)

    # Membrane parameters
    water_content = wt_percent_membrane*EW/M_H2O
    conductivity_membrane = (0.005139*water_content-0.00326)*exp(1268*(1/303-1/T_op_K))
    resistance_membrane = (t_membrane*1e-4)/conductivity_membrane     # units in ohm.cm2

    j0 = Loading_IrO2*Mass_activity_IrO2*exp(-j0_ref_E_act/R*(1/T_op_K-1/(St_T+Kelvin)))

    def current_density(eta):
        A = j0*(exp(((1-beta)*n*F*eta)/(R*(T_op_K)))-exp((-beta*n*F*eta)/(R*(T_op_K))))
        return A/(1+A/j_lim)

    V_t = list(linspace(Reversible_Voltage,V_final,100))

    def solve_eta(eta, V_t):
        A = j0*(exp(((1-beta)*n*F*eta)/(R*(T_op_K)))-exp((-beta*n*F*eta)/(R*(T_op_K))))
        j = A/(1+A/j_lim)
        eta_mt = ((R*T_op_K)*log(j_lim/(j_lim-j)))/(n*F)
        eta_ohm = j*resistance_membrane
        return Reversible_Voltage + eta + eta_ohm + eta_mt - V_t

    sol_eta = list(map(lambda vt:root_scalar(solve_eta,x0=0, x1=0.1, method='secant',args=(vt)).root,V_t))
    sol_current = list(map(lambda k:current_density(k),sol_eta))
    eta_mt = list(map(lambda y:((R*T_op_K)*log(j_lim/(j_lim-y)))/(n*F),sol_current))
    eta_ohm = list(map(lambda y:y*resistance_membrane,sol_current))

    return V_t,sol_current,sol_eta,eta_mt,eta_ohm,Reaction_del_H,Thermoneutral_Voltage,Reaction_del_G,Reversible_Voltage,conductivity_membrane

input_dictionary = {"Operating Temperature":T_op,"Operating Pressure":P_op,"Final Operating Voltage":V_final,"Limiting Current Density":j_lim,"Mass Activity of IrO₂":Mass_activity_IrO2,"IrO₂ Loading":Loading_IrO2,"IrO₂ Activation Energy at 25°C":j0_ref_E_act,"Asymmetric Factor (ß)":beta,"Membrane Weight%":wt_percent_membrane,"Membrane Equivalent Weight(EW)":EW,"Membrane Thickness":t_membrane}


results = list(Electrolysis_simulation(*list(input_dictionary.values())))

output_dictionary = {"Potential (V)":results[0],"Current Density (j)":results[1],"Activation Overpotential (ɳ)":results[2],"Mass Transport Losses (ɳₘₜ)":results[3],"Ohmic Losses (ɳₒₕₘ)":results[4],"Enthalpy of Reaction":results[5],"Thermoneutral Voltage":results[6],"Gibbs Free Energy of Reaction":results[7],"Reversible Voltage":results[8],"Membrane Conductivity":results[9]}

#st.write(Electrolysis_simulation(*list(input_dictionary.values()))[1])

st.divider()
st.header("Single Value Results")
# -------------------------
# ISOLATED SESSION STATE
# -------------------------
if "scalar_section" not in st.session_state:
    st.session_state.scalar_section = {
        "values": [{"id": str(uuid.uuid4()), "value": ""}]
    }

scalar_state = st.session_state.scalar_section

# -------------------------
# BASE RESULTS
# -------------------------
single_value_keys = [
    "Enthalpy of Reaction",
    "Thermoneutral Voltage",
    "Gibbs Free Energy of Reaction",
    "Reversible Voltage",
    "Membrane Conductivity"
]

base_results = {key: output_dictionary[key] for key in single_value_keys}

# -------------------------
# MULTI-FACTOR OPTION
# -------------------------
multi_scalar = st.checkbox("Run multi-factor experiment for single value results?")

multi_values_scalar = []

if multi_scalar:

    multi_factor_scalar = st.radio(
        "Select parameter to vary",
        ["Operating Temperature", "Operating Pressure", "Final Operating Voltage",
         "Limiting Current Density", "Mass Activity of IrO₂", "IrO₂ Loading",
         "IrO₂ Activation Energy at 25°C", "Asymmetric Factor (ß)",
         "Membrane Weight%", "Membrane Equivalent Weight(EW)", "Membrane Thickness"],
        horizontal=True
    )

    # Add value button
    if st.button("➕ Add Value (Results Section)"):
        scalar_state["values"].append({
        "id": str(uuid.uuid4()),
        "value": ""
    })

    # Input fields
    for item in scalar_state["values"]:
        vid = item["id"]
    
        col1, col2 = st.columns([5,1])
    
        with col1:
            item["value"] = st.text_input(
                "Value",
                value=item["value"],
                key=f"scalar_val_{vid}"
            )
    
        with col2:
            if st.button("❌", key=f"remove_scalar_{vid}"):
    
                scalar_state["values"] = [
                    v for v in scalar_state["values"] if v["id"] != vid
                ]
    
                # clean session state safely
                keys_to_delete = [k for k in st.session_state.keys() if vid in k]
                for k in keys_to_delete:
                    del st.session_state[k]
    
                st.rerun()

    multi_values_scalar = [item["value"] for item in scalar_state["values"]]

# -------------------------
# BUILD TABLE
# -------------------------
df_results = pd.DataFrame({
    "Parameter": single_value_keys,
    "Base Case": [base_results[k] for k in single_value_keys]
})

# -------------------------
# ADD MULTI-FACTOR COLUMNS
# -------------------------
if multi_scalar:

    base_input_dict = input_dictionary.copy()

    for val in multi_values_scalar:

        if val == "":
            continue

        val = float(val)

        if multi_factor_scalar == "Operating Pressure" and val == 0:
            val = 1

        # modify input
        temp_input = base_input_dict.copy()
        temp_input[multi_factor_scalar] = val

        results = list(Electrolysis_simulation(*list(temp_input.values())))
        temp_output = dict(zip(output_dictionary.keys(), results))

        col_name = f"{multi_factor_scalar}={val}"

        df_results[col_name] = [temp_output[k] for k in single_value_keys]

# -------------------------
# OPTIONAL POLISH
# -------------------------
df_results.iloc[:, 1:] = df_results.iloc[:, 1:].apply(
    lambda col: col.map(lambda x: round(x, 5))
)
units = {
    "Enthalpy of Reaction": "J/mol",
    "Thermoneutral Voltage": "V",
    "Gibbs Free Energy of Reaction": "J/mol",
    "Reversible Voltage": "V",
    "Membrane Conductivity": "S/cm"
}

df_results["Unit"] = df_results["Parameter"].map(units)
# -------------------------
# DISPLAY
# -------------------------
st.dataframe(df_results, use_container_width=True,hide_index=True)




st.divider()
st.header("Plot Configuration")


if "plots" not in st.session_state:
    st.session_state.plots = []

if "multi_values" not in st.session_state:
    st.session_state.multi_values = {}

# -------------------------
# ADD NEW PLOT BUTTON
# -------------------------
if st.button("Add Plot"):
    st.session_state.plots.append({
        "id": str(uuid.uuid4()),
        "x": None,
        "y": [],
        "multi_exp": False,
        "multi_factor": None,
        "values": []
    })
# -------------------------
# LOOP THROUGH PLOTS
# -------------------------
base_input_dict = input_dictionary.copy()  # preserve initial inputs

for idx, plot in enumerate(st.session_state.plots):
    pid = plot.get("id", str(uuid.uuid4()))
    col_title, col_remove = st.columns([6, 1])

    with col_title:
        st.subheader(f"Plot {idx+1}")

    with col_remove:
        remove_clicked = st.button(
            "❌",
            key=f"remove_{pid}",   # use unique ID
            disabled=len(st.session_state.plots) == 1
        )
        if remove_clicked:
            keys_to_delete = [k for k in st.session_state.keys() if pid in k]
            for k in keys_to_delete:
                del st.session_state[k]
            st.session_state.plots.pop(idx)
            st.rerun()

    fig, ax = subplots()

    # -------------------------
    # AXIS SELECTION
    # -------------------------
    col1, col2 = st.columns(2)
    with col1:
        plot["x"] = st.selectbox(
            "Choose X-axis",
            ["Potential (V)", "Current Density (j)", "Activation Overpotential (ɳ)",
             "Mass Transport Losses (ɳₘₜ)", "Ohmic Losses (ɳₒₕₘ)"],
            key=f"x_{pid}"
        )
    with col2:
        plot["y"] = st.multiselect(
            "Choose Y-axis",
            ["Potential (V)", "Current Density (j)", "Activation Overpotential (ɳ)",
             "Mass Transport Losses (ɳₘₜ)", "Ohmic Losses (ɳₒₕₘ)"],
            key=f"y_{pid}"
        )

    # -------------------------
    # MULTI EXPERIMENT
    # -------------------------
    plot["multi_exp"] = st.checkbox(
        "Run multi-factor experiment?",
        value=plot.get("multi_exp", False),
        key=f"multi_{pid}"
    )
    if plot["multi_exp"]:
        plot["multi_factor"] = st.radio(
            "Choose parameter",
            ["Operating Temperature", "Operating Pressure", "Final Operating Voltage",
             "Limiting Current Density", "Mass Activity of IrO₂", "IrO₂ Loading",
             "IrO₂ Activation Energy at 25°C", "Asymmetric Factor (ß)",
             "Membrane Weight%", "Membrane Equivalent Weight(EW)", "Membrane Thickness"],
            key=f"factor_{pid}", horizontal=True
        )

        col_add_val, _ = st.columns([5,1])
        with col_add_val:
            if st.button(f"➕ Add Value", key=f"add_val_{pid}"):
                plot.setdefault("values", []).append("")

        # Add remove buttons for each value
        for v_idx, val in enumerate(plot.get("values", [])):
            col_val, col_remove_val = st.columns([5, 1])
            with col_val:
                plot["values"][v_idx] = st.text_input(
                    f"Value {v_idx+1}",
                    value=plot["values"][v_idx],
                    key=f"value_{pid}_{v_idx}"
                )
            with col_remove_val:
                if st.button("❌", key=f"remove_val_{pid}_{v_idx}"):
                    plot["values"].pop(v_idx)
                    keys_to_delete = [k for k in st.session_state.keys() if f"value_{pid}_{v_idx}" in k]
                    for k in keys_to_delete:
                        del st.session_state[k]
                    st.rerun()

    # -------------------------
    # PLOTTING
    # -------------------------
    if plot["x"] and plot["y"]:
        x_data = output_dictionary[plot["x"]]

        # Set y-axis label: first y-parameter if multiple, else the single one
        ax.set_ylabel(plot["y"][0] if len(plot["y"]) > 1 else plot["y"][0])
        ax.set_xlabel(plot["x"])

        # --- Base case curves ---
        for y in plot["y"]:
            y_data = output_dictionary[y]
            # legend = parameter name if multiple y's, else "Base Case"
            label = y if len(plot["y"]) > 1 else "Base Case"
            ax.plot(x_data, y_data, label=label, linewidth=2)

        # --- Multi-factor experiment curves ---
        if plot.get("multi_exp", False):
            for val in plot.get("values", []):
                if val == "":
                    continue
                val = float(val)
                if plot["multi_factor"] == "Operating Pressure" and val == 0:
                    val = 1

                temp_input = base_input_dict.copy()
                temp_input[plot["multi_factor"]] = val
                results = list(Electrolysis_simulation(*list(temp_input.values())))
                temp_output = dict(zip(output_dictionary.keys(), results))

                for y in plot["y"]:
                    y_data = temp_output[y]
                    if len(plot["y"]) > 1:
                        label_name = f"{y} ({plot['multi_factor']}={val})"
                    else:
                        label_name = f"{plot['multi_factor']}={val}"
                    ax.plot(x_data, y_data, label=label_name)

            ax.legend()

    st.pyplot(fig)

    # -------------------------
    # ADD ANOTHER PLOT BUTTON AT THE END
    # -------------------------
    if idx == len(st.session_state.plots) - 1:
        if st.button("➕ Add Another Plot"):
            st.session_state.plots.append({
                "id": str(uuid.uuid4()),
                "x": None,
                "y": [],
                "multi_exp": False,
                "multi_factor": None,
                "values": []
            })
            st.rerun()
