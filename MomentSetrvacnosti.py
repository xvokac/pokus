from PyQt5.QtWidgets import QApplication, QComboBox, QTextEdit, QMessageBox, QLabel, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QTableWidget, QTableWidgetItem
from PyQt5.QtGui import QClipboard
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import math
import numpy as np
from datetime import datetime


# Funkce pro ArcTan2
def arctan2(y, x):
    pi = np.pi
    pi_2 = np.pi / 2

    if x > 0:
        return np.arctan(y / x)
    elif x < 0:
        return np.arctan(y / x) + pi * np.sign(y)
    elif x == 0:
        return pi_2 * np.sign(y)

#Funkce pro převod čísla do inenýrské notace
def to_engineering_notation(value):
    if value == 0:
        return "0.00e+00"
    
    exponent = int(math.floor(math.log10(abs(value)) // 3 * 3))  # Zajistí, že exponent je násobkem 3
    mantissa = value / (10 ** exponent)
    return f"{mantissa:.3f}e{exponent:+03d}"

#Výpočet plochy průřezu a těžiště
def cal_area_centrum(x,y):
    A = 0
    xT = 0
    yT = 0
    # Výpočet plochy a těžiště
    for i in range(len(x) - 1):
        A += (x[i] * y[i + 1] - x[i + 1] * y[i]) / 2
        xT += (x[i] + x[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i])
        yT += (y[i] + y[i + 1]) * (x[i] * y[i + 1] - x[i + 1] * y[i])
    xC = xT / (6 * A)
    yC = yT / (6 * A)
    return A, xC, yC

def cal_MoI(x,y):
    # Výpočet momentů setrvačnosti k ose x a y
    Ix = 0
    Iy = 0
    Ixy = 0
    for i in range(len(x) - 1):
        Ix += ((y[i]**2 + y[i] * y[i + 1] + y[i + 1]**2) * (x[i] * y[i + 1] - x[i + 1] * y[i])) / 12
        Iy += ((x[i]**2 + x[i] * x[i + 1] + x[i + 1]**2) * (x[i] * y[i + 1] - x[i + 1] * y[i])) / 12
        Ixy += -((y[i] - y[i + 1]) * (3 * x[i]**2 * y[i] + x[i]**2 * y[i + 1] + x[i + 1]**2 * y[i] + 
                                      3 * x[i + 1]**2 * y[i + 1] + 2 * x[i] * x[i + 1] * y[i] + 
                                      2 * x[i] * x[i + 1] * y[i + 1])) / 24
    return Ix, Iy, Ixy

def cal_cetr_MoI(Ix, Iy, Ixy):
    # Úhel hlavních os
    if Ixy == 0:
        alpha = 0
    else:
        alpha = 0.5 * arctan2(2 * Ixy, Iy - Ix)
    
    # Hlavní momenty setrvačnosti
    IxC = Ix * np.cos(alpha)**2 + Iy * np.sin(alpha)**2 - Ixy * np.sin(2 * alpha)
    IyC = Ix * np.sin(alpha)**2 + Iy * np.cos(alpha)**2 + Ixy * np.sin(2 * alpha)
    return IxC, IyC, alpha


class MainWindow(QWidget):
    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()

        # Nastavení názvu okna
        self.setWindowTitle("Výpočet momentů setrvačnosti plochy")

        # Layout pro výběr jednotky délky
        length_unit_layout = QHBoxLayout()

        # Přidání popisku pro jednotku délky
        self.unit_label = QLabel("Jednotka délky:")
        length_unit_layout.addWidget(self.unit_label)

        # Přidání comboboxu pro výběr jednotky délky
        self.unit_combobox = QComboBox()
        self.unit_combobox.addItems(["mm", "cm", "m", "in", "ft"])
        length_unit_layout.addWidget(self.unit_combobox)

        # Přidání layoutu s jednotkou délky do hlavního layoutu
        self.layout.addLayout(length_unit_layout)

        # Popisek nad tabulkou
        self.label = QLabel("Tabulka souřadnic (X, Y) uzavřeného polygonu orientovaného proti směru hodinových ručiček:")
        self.layout.addWidget(self.label)  # Přidání popisku nad tabulku

        # Tabulka pro souřadnice
        self.table = QTableWidget(0, 2)  # 0 řádků, 2 sloupce (X, Y)
        self.table.setHorizontalHeaderLabels(["X", "Y"])
        self.layout.addWidget(self.table)


        # Tlačítko pro přidání řádku
        self.add_row_button = QPushButton("Přidat řádek")
        self.add_row_button.clicked.connect(self.add_row)
        self.layout.addWidget(self.add_row_button)

        # Horizontální layout pro tlačítka
        control_layout = QHBoxLayout()  


        # Tlačítko pro odstranění vybraného řádku
        self.remove_row_button = QPushButton("Odstranit vybraný řádek")
        self.remove_row_button.clicked.connect(self.remove_selected_row)
        control_layout.addWidget(self.remove_row_button)

        # Tlačítko "Vymazat tabulku" pro provedení výpočtu
        self.clear_button = QPushButton("Vymazat tabulku")
        self.clear_button.clicked.connect(self.clear)
        control_layout.addWidget(self.clear_button)

        # Tlačítko pro vložení dat ze schránky
        self.paste_button = QPushButton("Paste from Clipboard")
        self.paste_button.clicked.connect(self.paste_from_clipboard)
        control_layout.addWidget(self.paste_button)

        # Přidání control_layout do hlavního layoutu
        self.layout.addLayout(control_layout)  

        # Tlačítko "Calculate" pro provedení výpočtu
        self.calculate_button = QPushButton("Výpočet")
        self.calculate_button.clicked.connect(self.calculate)
        self.layout.addWidget(self.calculate_button)

        # Přidání textového pole pro výstupy
        self.output_text = QTextEdit()
        self.output_text.setReadOnly(True)  # Nastavení jen pro čtení
        self.layout.addWidget(self.output_text)  # Přidání do hlavního layoutu

        # Nastavení hlavního layoutu
        self.setLayout(self.layout)

    def add_row(self):
        row_position = self.table.rowCount()
        self.table.insertRow(row_position)

        # Příklady pro přidání dat do buněk (můžete nechat prázdné)
        self.table.setItem(row_position, 0, QTableWidgetItem("0"))
        self.table.setItem(row_position, 1, QTableWidgetItem("0"))

    def remove_selected_row(self):
        """Odstraní vybraný řádek"""
        selected_row = self.table.currentRow()  # Získá aktuálně vybraný řádek
        if selected_row >= 0:
            self.table.removeRow(selected_row)

    def clear(self):
        # Vymazání tabulky před vložením nových dat
        self.table.setRowCount(0)

    def paste_from_clipboard(self):
        clipboard = QApplication.clipboard()
        data = clipboard.text()  # Získá textová data ze schránky (CSV formát)

        # Rozdělení dat podle řádků a jejich filtrace
        rows = [row for row in data.split('\n') if row.strip()]  # Ignoruje prázdné řádky

        for row_data in rows:
            columns = row_data.split('\t')  # Rozdělení na sloupce pomocí tabulátorů

            # Ověření správného počtu sloupců
            if len(columns) < 2:
                QMessageBox.warning(self, "Chyba", "Každý řádek musí obsahovat dvě hodnoty (X a Y).")
                continue  # Přeskočí řádky s nesprávným počtem hodnot

            # Přidání nového řádku do tabulky
            row_position = self.table.rowCount()
            self.table.insertRow(row_position)

            for column_position, cell_data in enumerate(columns[:2]):  # Pouze první dva sloupce
                # Nahradí čárku za tečku, pokud je to potřeba
                cell_data = cell_data.replace(",", ".").strip()
                if cell_data:  # Zajistí, že prázdné hodnoty nejsou vkládány
                    self.table.setItem(row_position, column_position, QTableWidgetItem(cell_data))
                else:
                    self.table.setItem(row_position, column_position, QTableWidgetItem("0"))  # Např. nahradí prázdné buňky nulou

    def append_output(self, text):
        """Funkce pro přidání textu do výstupního pole."""
        self.output_text.append(text)

    def get_selected_unit(self):
        """Funkce pro získání aktuálně vybrané jednotky."""
        return self.unit_combobox.currentText()




    def calculate(self):
        # Zavřít předchozí okno grafu, pokud existuje
        plt.close('all')
        # Vymazat obsah textového pole před novým výpočtem
        self.output_text.clear()

        """Kontroluje, zda jsou hodnoty v tabulce čísla, a provádí výpočet"""
        num_rows = self.table.rowCount()

        x_values = []
        y_values = []
        
        if num_rows < 3:
            # Pokud je málo bodů, zobrazit varování
            QMessageBox.warning(self, "Chyba", "Minimální počet řádků je 3.")
            return

        # Procházení všech buněk tabulky
        for row in range(num_rows):
            try:
                x_item = self.table.item(row, 0).text()
                y_item = self.table.item(row, 1).text()

                # Zkontrolovat, jestli jsou položky čísla
                x_value = float(x_item)
                y_value = float(y_item)

                x_values.append(x_value)
                y_values.append(y_value)

            except ValueError:
                # Pokud není číslo, zobrazit varování
                QMessageBox.warning(self, "Chyba", f"Řádek {row + 1}: Hodnoty musí být čísla.")
                return
        
        # Uzavrit polygon - pokud není v tabulce, tak se uzavře
        if x_values[0] != x_values[num_rows-1] or y_values[0] != y_values[num_rows-1]:
            x_values = np.append(x_values, x_values[0])
            y_values = np.append(y_values, y_values[0])
       

        """
        Vlastní výpočet momentů setrvačnosti
        """
        # Plocha průřezu a souřadnice těžiště
        A, xC, yC = cal_area_centrum(x_values,y_values)
        # Momenty setrvačnosti k osám x a y, v kterých jsou zadané souřadnice
        Ix, Iy, Ixy = cal_MoI(x_values,y_values)

        # Výpočet těžišťových momentů setrvačnosti - osy rovnoběžné s původním souř. sys.
        Ix = Ix - A * yC**2  # Steinerova věta
        Iy = Iy - A * xC**2
        Ixy = Ixy - A * xC * yC

        # Výpočet hlavních centrálních mometů setrvačnosti
        IxC, IyC, alpha = cal_cetr_MoI(Ix, Iy, Ixy)

        # Poloměry setrvačnosti
        ixc = np.sqrt(IxC/A)
        iyc = np.sqrt(IyC/A)

        """
        ############################################################
        # Vykreslení  v grafu
        """

        selected_unit = self.get_selected_unit()


        # Délka šipek hl. os souřadnic - odhad optimalni velikosti
        arrow_length = np.max([np.abs(np.max(x_values)-np.mean(x_values)), np.abs(np.min(x_values)-np.mean(x_values)),
                       np.abs(np.max(y_values)-np.mean(y_values)), np.abs(np.min(y_values)-np.mean(y_values)) ]) 

        # Složky šipek pro osy x a y
        arrow_x_dx = arrow_length * np.cos(alpha)  # Složka ve směru osy x
        arrow_x_dy = arrow_length * np.sin(alpha)  # Složka ve směru osy y

        arrow_y_dx = -arrow_length * np.sin(alpha)  # Složka pro kolmicu ve směru y
        arrow_y_dy = arrow_length * np.cos(alpha)   # Složka pro kolmicu ve směru y
        
        # Parametry elipsy setrvačnosti
        # Generování bodů elipsy
        t = np.linspace(0, 2 * np.pi, 100)  # Parametrický úhel
        
        # Parametrické vyjádření elipsy (v počátku)
        x_ellipse = iyc * np.cos(t)
        y_ellipse = ixc * np.sin(t)

        # Natočení elipsy o úhel alpha_rad
        x_rot = xC + (x_ellipse * np.cos(alpha) - y_ellipse * np.sin(alpha))
        y_rot = yC + (x_ellipse * np.sin(alpha) + y_ellipse * np.cos(alpha))
       
        # Vykreslení grafu
        plt.plot(x_values, y_values, '-o', label='Uzavřený polygon')  # Linie s body
        plt.fill(x_values, y_values, 'b', alpha=0.2)  # Vyplnění uzavřené oblasti
        
        # Vykreslení souřadnicových os (šipek) v bodě (xC, yC)
        plt.quiver(xC, yC, arrow_x_dx, arrow_x_dy, angles='xy', scale_units='xy', scale=1, color='r', label="Osa $x_c$")  # Osa X
        plt.quiver(xC, yC, arrow_y_dx, arrow_y_dy, angles='xy', scale_units='xy', scale=1, color='g', label="Osa $y_c$")  # Osa Y
        
        # Vykreslení elipsy
        plt.plot(x_rot, y_rot, label="Elipsa setrvačnosti", color='b')
        
        # Vykreslení středu elipsy
        plt.plot(xC, yC, 'bo', label="Těžiště")
        
        # Nastavení grafu
        plt.title("Průřez")
        plt.xlabel(f"$X$ [{selected_unit}]")
        plt.ylabel(f"$Y$ [{selected_unit}]")
        plt.legend()
        plt.grid(True)
        plt.axis('equal')  # Zachování poměru os pro správné vykreslení

        plt.show()

        """
        #Výpis hodnot #################################################################################
        """
        now = datetime.now() # Získání aktuálního data a času
        date_time_str = now.strftime("%d-%m-%Y %H:%M:%S")
        self.append_output(f'##### Výpočet plošných momentů setrvačnosti ({date_time_str}) #####')
        self.append_output('# Průřez zadán uzavřeným polygonem:')
        self.append_output(f'X [{selected_unit}]\tY [{selected_unit}]')
        for row in range(num_rows):
            self.append_output(f'{x_values[row]}\t{y_values[row]}')
        self.append_output('# Výsledky výpočtu:')
        self.append_output(f'Pocha průřezu\t\t\tA =\t{to_engineering_notation(A)}\t{selected_unit}^2')
        self.append_output(f'Souřadnice těžiště\t\tX_t =\t{to_engineering_notation(xC)}\t{selected_unit}')
        self.append_output(f'Souřadnice těžiště\t\tY_t =\t{to_engineering_notation(yC)}\t{selected_unit}')
        self.append_output(f'Moment setrvačnosti k těžišťové ose\tI_x =\t{to_engineering_notation(Ix)}\t{selected_unit}^4')
        self.append_output(f'Moment setrvačnosti k těžišťové ose\tI_y =\t{to_engineering_notation(Iy)}\t{selected_unit}^4')
        self.append_output(f'Deviační moment k těžišťovým osám\tI_xy =\t{to_engineering_notation(Ixy)}\t{selected_unit}^4')
        self.append_output(f'Natočení hlavních os setrvačnosti\talpha_c =\t{to_engineering_notation(alpha*180/np.pi)}\tdeg')
        self.append_output(f'Hlavní centrální moment setrvačnosti\tI_xc =\t{to_engineering_notation(IxC)}\t{selected_unit}^4')
        self.append_output(f'Hlavní centrální moment setrvačnosti\tI_yc =\t{to_engineering_notation(IyC)}\t{selected_unit}^4')
        self.append_output(f'Poloměr setrvačnosti průřezu\t\ti_xc =\t{to_engineering_notation(ixc)}\t{selected_unit}')
        self.append_output(f'Poloměr setrvačnosti průřezu\t\ti_yc =\t{to_engineering_notation(iyc)}\t{selected_unit}')
        self.append_output(' ')


        
# Aplikace
app = QApplication([])
window = MainWindow()
window.show()
sys.exit(app.exec_())
