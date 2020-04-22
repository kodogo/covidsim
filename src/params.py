N = 5000000    # Population
X = 0    # Initial infections
psi = 10    # New infected travellers per day
t_max = 0    # Time of highest transmission potential
q_max = 1500    # Maximum isolation capacity
t_iso1 = 0    # Isolation start time
iso_months = 0    # Isolation months
t_iso2 = t_iso1 + (iso_months * 365.25 / 12)    # Isolation end time
c_home = 0.5    # Contact reduction for home isolation
c_cont = 0.25    # General contact reduction
t_cont1 = 0    # Contact reduction start time
contact_months = 0    # Contact reduction months
t_cont2 = t_cont1 + (contact_months * 365.25 / 12)    # Contact reduction end time
R0 = 3.5    # Basic reproduction number
a = 0.25    # Amplitude of seasonal variation
a_max = 105    # Peak date of seasonal variation - 15/7/2020 for 1 April start date
D_e = 4    # Duration of latent period
n_e = 16    # Stages for latent period
epsilon = n_e / D_e    # Stage rate for latent period
D_p = 1    # Duration of prodromal period
n_p = 16    # Stages for prodromal period
phi = n_p / D_p    # Stage rate for prodromal period
i_p = 0.5    # Relative infectiosness in prodromal period
D_i = 10    # Duration of symptomatic period
n_i = 16    # Stages for infectious period
gamma = n_i / D_i    # Stage rate for infectios period
beta = R0/(i_p*D_p + D_i)   # Base contact rate
p_sick = 0.67    # Fraction of infected who become sick
p_consult = 0.4    # Fraction of infected who seek medical help
p_hosp = 0.044    # Fraction of infected who are hospitalised
p_icu = 0.25    # Fraction of infected who are admitted to ICU
p_death = 0.0083    # Fraction of infected who die