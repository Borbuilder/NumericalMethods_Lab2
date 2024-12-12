import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.ticker as ticker
import numpy as np
from Lab_2_Module import Solver, MODE


class Window:
    def __init__(self, master):
        # variables
        self.solve_mode = MODE.TEST
        self.nodes_num = None
        self.solver = None
        self.table_data = None
        self.graph_mode = 0

        # window master
        self.master = master
        self.master.title("Лабораторная работа №2")
        self.master.geometry("1300x600")

        # Виджет выбора типа задачи
        self.mode_label = tk.Label(self.master, text="Выберите тип задачи:")
        self.mode_label.grid(row=0, column=0, sticky='w')

        # Dropdown for mode selection
        self.mode_selector = ttk.Combobox(self.master, values=["Тестовая", "Основная"])
        self.mode_selector.grid(row=1, column=0, sticky='we')
        self.mode_selector.set("Тестовая")
        self.mode_selector.bind("<<ComboboxSelected>>", self.on_mode_selected)

        # Display widgets for "Тестовая задача" by default
        self.create_test_mode_widgets()

    def on_mode_selected(self, event):
        selected_mode = self.mode_selector.get()
        self.clear_widgets()

        # Choose which widgets to display based on selected task
        if selected_mode == "Тестовая":
            self.solve_mode = MODE.TEST
            self.create_test_mode_widgets()
        elif selected_mode == "Основная":
            self.solve_mode = MODE.MAIN
            self.create_main_mode_widgets()

    def clear_widgets(self):
        # Remove all widgets except the dropdown
        for widget in self.master.winfo_children():
            if widget not in [self.mode_selector, self.mode_label]:
                widget.destroy()

        # Clear variables
        self.nodes_num = None
        self.solver = None
        self.table_data = None
        self.graph_mode = 0

    def create_test_mode_widgets(self):
        self.master.grid_columnconfigure(0, weight=100)
        self.master.grid_columnconfigure(1, weight=2)
        self.master.grid_columnconfigure(2, weight=2)
        self.master.grid_columnconfigure(3, weight=2)

        # Input fields for test task
        tk.Label(self.master, text="Кол-во узлов:").grid(row=0, column=1, sticky='e')

        # Input fields
        self.entry_nodes_num = tk.Entry(self.master)

        # Set default values
        self.entry_nodes_num.insert(0, "1000")

        # Placement of input fields
        self.entry_nodes_num.grid(row=0, column=2, sticky='ew')

        # Calculate button
        self.calc_button = tk.Button(self.master, text="Вычислить", command=self.calculate)
        self.calc_button.grid(row=6, column=0, sticky='we')

        # Create Treeview to display the table
        self.table = ttk.Treeview(self.master, columns=["Col" + str(i) for i in range(5)], show="headings")
        self.table.grid(row=7, rowspan=2, column=0, columnspan=3, sticky='nsew', padx=5, pady=5)

        # Add horizontal and vertical scrolling
        vsb = ttk.Scrollbar(self.master, orient="vertical", command=self.table.yview)
        vsb.grid(row=7, column=2, sticky='sne', rowspan=3)
        self.table.configure(yscrollcommand=vsb.set)

        hsb = ttk.Scrollbar(self.master, orient="horizontal", command=self.table.xview)
        hsb.grid(row=9, column=0, columnspan=3, sticky='esw')
        self.table.configure(xscrollcommand=hsb.set)

        # Configure stretching
        self.master.grid_rowconfigure(8, weight=1)  # Allows the table to expand vertically
        self.master.grid_columnconfigure(0, weight=1)  # Allows the table to expand horizontally

        # Create figure for the plot
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.master)
        self.canvas.get_tk_widget().grid(row=6, column=4, rowspan=5, sticky='se', padx=5, pady=5)

        # Create a Text widget for final report
        self.final_report_text = tk.Text(self.master, height=10, width=50)
        self.final_report_text.grid(row=0, column=4, rowspan=7, sticky='nwe', padx=5, pady=5)

        # Button for plotting U-V(x)
        self.plot_button = tk.Button(self.master, text="Построить график разности U(x) и V(x)", command=self.toggle_plot)
        self.plot_button.grid(row=7, column=4, sticky='we', padx=5, pady=5)

    def create_main_mode_widgets(self):
        self.master.grid_columnconfigure(0, weight=100)
        self.master.grid_columnconfigure(1, weight=2)
        self.master.grid_columnconfigure(2, weight=2)
        self.master.grid_columnconfigure(3, weight=2)

        # Input fields for test task
        tk.Label(self.master, text="Кол-во узлов:").grid(row=0, column=1, sticky='e')

        # Input fields
        self.entry_nodes_num = tk.Entry(self.master)

        # Set default values
        self.entry_nodes_num.insert(0, "1000")

        # Placement of input fields
        self.entry_nodes_num.grid(row=0, column=2, sticky='ew')

        # Calculate button
        self.calc_button = tk.Button(self.master, text="Вычислить", command=self.calculate)
        self.calc_button.grid(row=6, column=0, sticky='we')

        # Create Treeview to display the table
        self.table = ttk.Treeview(self.master, columns=["Col" + str(i) for i in range(5)], show="headings")
        self.table.grid(row=7, rowspan=2, column=0, columnspan=3, sticky='nsew', padx=5, pady=5)

        # Add horizontal and vertical scrolling
        vsb = ttk.Scrollbar(self.master, orient="vertical", command=self.table.yview)
        vsb.grid(row=7, column=2, sticky='sne', rowspan=3)
        self.table.configure(yscrollcommand=vsb.set)

        hsb = ttk.Scrollbar(self.master, orient="horizontal", command=self.table.xview)
        hsb.grid(row=9, column=0, columnspan=3, sticky='esw')
        self.table.configure(xscrollcommand=hsb.set)

        # Configure stretching
        self.master.grid_rowconfigure(8, weight=1)  # Allows the table to expand vertically
        self.master.grid_columnconfigure(0, weight=1)  # Allows the table to expand horizontally

        # Create figure for the plot
        self.figure = Figure(figsize=(6, 4), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.master)
        self.canvas.get_tk_widget().grid(row=6, column=4, rowspan=5, sticky='se', padx=5, pady=5)

        # Create a Text widget for final report
        self.final_report_text = tk.Text(self.master, height=10, width=50)
        self.final_report_text.grid(row=0, column=4, rowspan=7, sticky='nwe', padx=5, pady=5)

        # Button for plotting U-V(x)
        self.plot_button = tk.Button(self.master, text="Построить график разности V(x) и V2(x)", command=self.toggle_plot)
        self.plot_button.grid(row=7, column=4, sticky='we', padx=5, pady=5)

    def calculate(self):
        try:
            # Get values from input fields
            self.nodes_num = int(self.entry_nodes_num.get())

            # Checks for correctness of values
            if self.nodes_num <= 1:
                messagebox.showerror("Ошибка", "Количество узлов должно быть больше 1.")
                return

            # Solving based on selected item in ComboBox
            self.solver = Solver()
            self.solver.Solve(self.nodes_num, self.solve_mode)
            self.table_data = self.solver.get_table()
            self.update_table()

            # Reference
            self.show_final_reference()

            # Graph
            if self.graph_mode == 0:
                if self.solve_mode == MODE.TEST:
                    self.plot_graph_UV_X()
                elif self.solve_mode == MODE.MAIN:
                    self.plot_graph_VV2_X()

            elif self.graph_mode == 1:
                if self.solve_mode == MODE.TEST:
                    self.plot_graph_UV_diference()
                elif self.solve_mode == MODE.MAIN:
                    self.plot_graph_VV2_diference()

        except ValueError:
            messagebox.showerror("Ошибка", "Пожалуйста, введите корректные значения.")

    def update_table(self):
        # Clear Treeview before updating
        for row in self.table.get_children():
            self.table.delete(row)

        # Set column headers
        if self.solve_mode == MODE.TEST:
            # columns = ["i", "X", "Num. sol.", "True sol.", "Difference"]
            columns = ["i", "xi", "U(xi)", "V(xi)", "U(xi) - V(xi)"]
            arr_width = [50, 100, 100, 100, 100]

        elif self.solve_mode == MODE.MAIN:
            columns = ["i", "xi", "V(xi)", "V2(xi)", "V(xi) - V2(xi)"]
            arr_width = [50, 100, 100, 100, 100]

        if list(self.table["columns"]) != columns:
            self.table["columns"] = columns
            for i in range(len(columns)):
                self.table.heading(columns[i], text=columns[i])
                self.table.column(columns[i], width=arr_width[i], minwidth=arr_width[i])  # Set column width

        # Add data to the table with formatting
        for row in self.table_data:
            formatted_row = [f"{value:.6g}" if isinstance(value, float) else value for value in row]
            self.table.insert("", "end", values=formatted_row)

    # Function for switching between graphs
    def toggle_plot(self):
        if self.graph_mode == 0:
            self.graph_mode = 1
            if self.solve_mode == MODE.TEST:
                self.plot_graph_UV_diference()
                self.plot_button.config(text="Построить графики U(x) и V(x)")
            elif self.solve_mode == MODE.MAIN:
                self.plot_graph_VV2_diference()
                self.plot_button.config(text="Построить график V(x) и V2(x)")

        elif self.graph_mode == 1:
            self.graph_mode = 0
            if self.solve_mode == MODE.TEST:
                self.plot_graph_UV_X()
                self.plot_button.config(text="Построить график разности U(x) и V(x)")
            elif self.solve_mode == MODE.MAIN:
                self.plot_graph_VV2_X()
                self.plot_button.config(text="Построить график разности V(x) и V2(x)")

    def plot_graph_UV_X(self):
        # Clearing the graph
        self.ax.clear()

        if self.solve_mode == MODE.TEST:
            # Extract X, U and V values from table data for the plot
            X_values = [row[1] for row in self.table_data]
            U_values = [row[2] for row in self.table_data]
            V_values = [row[3] for row in self.table_data]

            # Plot U(x) and V(x)
            self.ax.plot(X_values, U_values, label=f'U(x)', color='blue', alpha=0.7)
            self.ax.plot(X_values, V_values, label=f'V(x)', color='red', alpha=0.7)

        # Customize the plot
        self.ax.set_xlabel('X')
        # self.ax.xaxis.set_tick_params(labelsize=8)
        # self.ax.yaxis.set_tick_params(labelsize=8)
        self.ax.legend()
        self.ax.grid()

        # Автоматическая настройка размещения элементов графика
        self.ax.figure.tight_layout()

        # Update the plot
        self.canvas.draw()

    def plot_graph_UV_diference(self):
        # Clearing the graph
        self.ax.clear()

        if self.solve_mode == MODE.TEST:
            # Extract X, U and V values from table data for the plot
            X_values = [row[1] for row in self.table_data]
            UV_dif_values = [row[4] for row in self.table_data]

            # Plot U(x) and V(x)
            self.ax.plot(X_values, UV_dif_values, label=f'U(x)-V(x)', color='blue', alpha=0.7)

        # Customize the plot
        self.ax.set_xlabel('X')
        # self.ax.xaxis.set_tick_params(labelsize=8)
        # self.ax.yaxis.set_tick_params(labelsize=8)
        self.ax.legend()
        self.ax.grid()

        # Автоматическая настройка размещения элементов графика
        self.ax.figure.tight_layout()

        # Update the plot
        self.canvas.draw()

    def plot_graph_VV2_X(self):
        # Clearing the graph
        self.ax.clear()

        if self.solve_mode == MODE.MAIN:
            # Extract X, U and V values from table data for the plot
            X_values = [row[1] for row in self.table_data]
            V_values = [row[2] for row in self.table_data]
            V2_values = [row[3] for row in self.table_data]

            # Plot U(x) and V(x)
            self.ax.plot(X_values, V_values, label=f'V(x)', color='blue', alpha=0.7)
            self.ax.plot(X_values, V2_values, label=f'V2(x)', color='red', alpha=0.7)

        # Customize the plot
        self.ax.set_xlabel('X')
        # self.ax.xaxis.set_tick_params(labelsize=8)
        # self.ax.yaxis.set_tick_params(labelsize=8)
        self.ax.legend()
        self.ax.grid()

        # Автоматическая настройка размещения элементов графика
        self.ax.figure.tight_layout()

        # Update the plot
        self.canvas.draw()

    def plot_graph_VV2_diference(self):
        # Clearing the graph
        self.ax.clear()

        if self.solve_mode == MODE.MAIN:
            # Extract X, U and V values from table data for the plot
            X_values = [row[1] for row in self.table_data]
            UV_dif_values = [row[4] for row in self.table_data]

            # Plot U(x) and V(x)
            self.ax.plot(X_values, UV_dif_values, label=f'V(x)-V2(x)', color='blue', alpha=0.7)

        # Customize the plot
        self.ax.set_xlabel('X')
        # self.ax.xaxis.set_tick_params(labelsize=8)
        # self.ax.yaxis.set_tick_params(labelsize=8)
        self.ax.legend()
        self.ax.grid()

        # Автоматическая настройка размещения элементов графика
        self.ax.figure.tight_layout()

        # Update the plot
        self.canvas.draw()

    def show_final_reference(self):
        # Clear the previous report text
        self.final_report_text.delete(1.0, tk.END)

        # Insert new information into the Text widget
        if self.solve_mode == MODE.TEST:
            info = (
                f"Для решения задачи использована равномерная сетка с числом\nразбиений n={self.nodes_num}.\n"
                f"Разбиение происходит на отрезке [{self.solver.get_diaposon()[0]}; {self.solver.get_diaposon()[1]}] с шагом {self.solver.get_step()}.\n"
                f"Задача должна быть решена с погрешностью не более \u03B5=0.5*10-6.\n"
                f"Задача решена с погрешностью \u03B51={self.solver.get_max_diff()}.\n"
                f"Максимальное отклонение аналитического и численного решений наблюдается в точке x={self.solver.get_x_diff()}.\n"
                f"Координаты точки разрыва x={self.solver.get_break_point()}."
            )
        else:
            info = (
                f"Для решения задачи использована равномерная сетка с числом\nразбиений n={self.nodes_num}.\n"
                f"Разбиение происходит на отрезке [{self.solver.get_diaposon()[0]}; {self.solver.get_diaposon()[1]}] с шагом {self.solver.get_step()}.\n"
                f"Задача должна быть решена с заданной точностью \u03B5=0.5*10-6.\n"
                f"Задача решена с точностью \u03B52={self.solver.get_max_diff()}.\n"
                f"Максимальная разность численных решений в общих узлах сетки наблюдается в точке x={self.solver.get_x_diff()}.\n"
                f"Координаты точки разрыва x={self.solver.get_break_point()}."
            )

        self.final_report_text.insert(tk.END, info)


def create_gui():
    root = tk.Tk()
    app = Window(root)
    root.mainloop()


def print_table(table):
    for i in table:
        for j in i:
            # print(str(j).ljust(20), end="\t")
            print(j, end="\t")

        print()


if __name__ == "__main__":
    create_gui()
