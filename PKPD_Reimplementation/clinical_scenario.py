def generate_infusion_func(weight_kg, bolus_mg_kg=2.0, bolus_time_min=1.0, maint_mg_kg_h=10.0, maint_time_min=60.0):
    """
    Returns a clinical infusion scenario schedule as a mathematically callable function.
    """
    bolus_dose_mg = bolus_mg_kg * weight_kg
    maint_rate_mg_min = (maint_mg_kg_h * weight_kg) / 60.0

    def infusion_rate(t):
        if t <= bolus_time_min:
            # Distribute bolus evenly over bolus_time_min
            return bolus_dose_mg / bolus_time_min
        elif t <= (bolus_time_min + maint_time_min):
            return maint_rate_mg_min
        else:
            return 0.0
            
    return infusion_rate
