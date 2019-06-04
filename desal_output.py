import desal_model


out1 = {
    "Na":   desal_model.output_water[0],
    "Ca":   desal_model.output_water[1],
    "Mg":   desal_model.output_water[2],
    "K":    desal_model.output_water[3],
    "Cl":   desal_model.output_water[4],
    "SO4":  desal_model.output_water[5],
    "CO3":  desal_model.output_water[6],
    "HCO3": desal_model.output_water[7],
    "NO3":  desal_model.output_water[8],
    "TDS":  sum(desal_model.output_water)
}

out1 = {x: round(out1[x], 1) for x in out1}
    
out2 = round(desal_model.tmp, 2)

out3 = {x: round(desal_model.concentrate1.si[x], 3) for x in desal_model.concentrate1.si}

out4 = round(desal_model.sec, 3)

out5 = round(desal_model.flux*3600, 1)

print("\nNumber of pass:")
print(desal_model.no_of_pass)

print("\nPermeate flux:")
print(out5, "LMH")

print("\nPermeate concentration:")
for x in out1:
    print(x, out1[x], "mg/L")

print("\nTransmembrane pressure:", out2, "bar")

print("\nSpecific energy consumption:", out4, "kWs/m3")

print("\nConcentrate SI:")
for x in out3:
    if out3[x] > 0:
        print(x, out3[x])
print("LSI", desal_model.concentrate1.scale[0])
print("SDSI", desal_model.concentrate1.scale[1])
print("\n")