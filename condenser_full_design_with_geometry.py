
import math
import streamlit as st
from CoolProp.CoolProp import PropsSI
from scipy.optimize import fsolve

st.title("Freon Condenser Design: Full Calculation with Coil Geometry")

st.sidebar.header("1. Refrigerant and Conditions")
fluid = st.sidebar.selectbox("Refrigerant", ["R134a", "R407C"])
m_dot_freon = st.sidebar.number_input("Freon Mass Flow Rate (kg/s)", value=0.6)
T_super = st.sidebar.number_input("Superheated Inlet Temp (°C)", value=95.0) + 273.15
T_cond = st.sidebar.number_input("Condensing Temp (°C)", value=57.0) + 273.15
T_sub = st.sidebar.number_input("Subcooled Outlet Temp (°C)", value=52.0) + 273.15
P_cond = st.sidebar.number_input("Condensing Pressure (Pa)", value=2352000.0)

st.sidebar.header("2. Air Side Conditions")
airflow = st.sidebar.number_input("Air Flow Rate (m³/s)", value=10.0)
air_density = 1.06
cp_air = 1006
air_mass_flow = airflow * air_density
face_velocity = st.sidebar.number_input("Air Face Velocity (m/s)", value=2.5)

T_air_sub = st.sidebar.number_input("Air Temp into Subcooling Zone (°C)", value=40.0) + 273.15
T_air_cond = st.sidebar.number_input("Air Temp into Condensing Zone (°C)", value=45.0) + 273.15
T_air_super = st.sidebar.number_input("Air Temp into Desuperheating Zone (°C)", value=52.0) + 273.15

st.sidebar.header("3. Heat Transfer U Values (W/m²·K)")
U_sub = st.sidebar.number_input("U for Subcooling", value=65.0)
U_cond = st.sidebar.number_input("U for Condensing", value=80.0)
U_desuper = st.sidebar.number_input("U for Desuperheating", value=70.0)

st.sidebar.header("4. Coil Geometry")
D_o = st.sidebar.number_input("Tube Outer Diameter (mm)", value=9.52) / 1000  # in meters
fpi = st.sidebar.number_input("Fins Per Inch (FPI)", value=10.0)
rows_max = st.sidebar.slider("Max Rows per Zone", min_value=2, max_value=6, value=4)

# Enthalpies
h_super = PropsSI("H", "T", T_super, "P", P_cond, fluid)
h_vap = PropsSI("H", "T", T_cond, "Q", 1, fluid)
h_liq = PropsSI("H", "T", T_cond, "Q", 0, fluid)
h_sub = PropsSI("H", "T", T_sub, "P", P_cond, fluid)

Q_desuper = m_dot_freon * (h_super - h_vap)
Q_cond = m_dot_freon * (h_vap - h_liq)
Q_sub = m_dot_freon * (h_liq - h_sub)

def solve_NTU_eps(eps_target, Cr):
    def eq(NTU):
        return 1 - math.exp((1 / Cr) * (NTU**0.22) * (math.exp(-Cr * NTU**0.78) - 1)) - eps_target
    return fsolve(eq, 1.0)[0]

def compute_area(Q, T_hot, T_air, m_dot_air, U, phase_change=False):
    C_air = m_dot_air * cp_air
    deltaT = T_hot - T_air
    eps = Q / (C_air * deltaT)
    NTU = -math.log(1 - eps) if phase_change else solve_NTU_eps(eps, Cr=1e-6)
    A = NTU * C_air / U
    T_air_out = T_air + Q / C_air
    return A, T_air_out - 273.15, NTU, eps

st.header("Zone-wise Heat Transfer and Geometry")

zones = [
    ("Subcooling", Q_sub, T_cond, T_air_sub, air_mass_flow, U_sub, False),
    ("Condensation", Q_cond, T_cond, T_air_cond, air_mass_flow, U_cond, True),
    ("Desuperheating", Q_desuper, T_super, T_air_super, air_mass_flow, U_desuper, False),
]

for name, Q, T_hot, T_air, m_dot_air, U, phase_change in zones:
    A, T_air_out, NTU, eps = compute_area(Q, T_hot, T_air, m_dot_air, U, phase_change)

    # Face area from airflow and face velocity
    face_area = airflow / face_velocity
    coil_width = math.sqrt(face_area)
    coil_height = face_area / coil_width

    fin_pitch = 1 / (fpi * 39.37)  # m between fins
    fin_rows = min(rows_max, int(A / (coil_width * fin_pitch)))
    tube_passes = max(1, int(fin_rows / 2))

    st.subheader(name)
    st.write(f"Heat Load: {Q/1000:.2f} kW")
    st.write(f"Required Area: {A:.2f} m²")
    st.write(f"NTU: {NTU:.2f}, Effectiveness: {eps:.3f}")
    st.write(f"Air Outlet Temp: {T_air_out:.2f} °C")
    st.write(f"Coil Face Area: {face_area:.3f} m² ({coil_width:.2f} × {coil_height:.2f} m)")
    st.write(f"Estimated Rows: {fin_rows}, Estimated Passes: {tube_passes}")
