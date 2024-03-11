Temp = ["1.00", "1.20", "1.40", "1.60", "1.80", "2.00", "2.20", "2.40", "2.60", "2.80", "3.00", "3.20", "3.40", "3.60", "3.80", "4.00", "4.20", "4.60", "4.80", "5.00"]

data = []
for i in range(len(Temp)):
    file = open(f"results-x/results-H0.01-T{Temp[i]}/data/gamma.dat", "r")
    lines = file.readlines()
    lines = lines[1:]
    file.close()
    for line in lines:
        numbers = line.split()
        list_numbers = []
        for number in numbers:
            list_numbers.append(float(number))
        data.append(list_numbers)

file = open("results/data/gamma.dat", "w")
file.write("#TEMP    Γ       p       γ                Δγ\n")
for i in range(len(data)):
    file.write(f"{data[i][0]:.2f}    {data[i][1]:.2f}    {data[i][2]:.2f}     {data[i][3]:.10f}     {data[i][4]:.10f}")
    file.write("\n")
file.close()

data = []
for i in range(len(Temp)):
    file = open(f"results-x/results-H0.01-T{Temp[i]}/data/mz.dat", "r")
    lines = file.readlines()
    lines = lines[1:]
    file.close()
    for line in lines:
        numbers = line.split()
        list_numbers = []
        for number in numbers:
            list_numbers.append(float(number))
        data.append(list_numbers)

file = open("results/data/mz.dat", "w")
file.write("#TEMP    Γ       p       <|mz|>\n")
for i in range(len(data)):
    file.write(f"{data[i][0]:.2f}    {data[i][1]:.2f}    {data[i][2]:.2f}     {data[i][3]:.10f}")
    file.write("\n")
file.close()