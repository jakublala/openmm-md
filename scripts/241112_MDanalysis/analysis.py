from src.analysis.system import plot_system_data

def main():
    date = "241111"
    project = "241109_INFconstruct"
    for protease in ["Q7", "Z1"]:
        for design in ["B30L4", "B30L7", "B30L10", "B40L10", "B40L10W", "B50W"]:
            filepath = f"../../data/{project}/output/{protease}-{design}/{date}/{protease}-{design}"
            plot_system_data(filepath)

if __name__ == "__main__":
    main()