################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 12.02.2021
#
#   Info: Function to show how Quantization works.
#   RoadMap: Use PyQT to create a GUI interface
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np


def plot_signals(freq_cont=0, freq_disc=0, quantization=0, points_per_cycle=5001):
    if freq_cont == 0 or freq_disc == 0 or quantization == 0:
        print("\t Erro! A frequencia de aquisição ou a quantização não foi definida!")
    else:
        t_cont = np.linspace(0.0, 1, num=points_per_cycle)
        y_cont = 2.5 * np.sin(2 * np.pi * freq_cont * t_cont) + 2.5
        dv = 5 / (2 ** quantization - 1)
        t_amos = linspace(0.0, 1, num=freq_disc + 1)
        y_amos = np.round((2.5*np.sin(2 * np.pi * freq_cont * t_amos) + 2.5)/dv) * dv
        # plots
        plt.step(t_amos, y_amos, 'o:', color=(0, 0.6, 0, 1), where='post', linewidth=2.0, markersize=12)
        plt.plot(t_cont, y_cont, 'b-', linewidth=3.0)
        plt.xlabel('Tempo [s]')
        plt.ylabel('Tensão [V]')
        if quantization <= 5:
            for i in range(2 ** quantization):
                plt.plot([0, 1], [i * dv,  i*dv], color=(0.7, 0.7, 0.7, 0.3), linewidth=1.0)

        if freq_disc <= 20:
            for t in t_amos:
                plt.plot([t, t], [0, 5], color=(0.7, 0.7, 0.7, 0.3), linewidth=1.0)
        # if (qt <= 50)
        #     for i=1:qt
        #     line([(i - 1) * tq(2)(i - 1) * tq(2)], [0 2], 'LineStyle', ':', 'Color', 0.85 * [1 1 1]);
        #     end
        # end

        plt.show()


def print_menu():
    print("\t As Opções são: ")
    print("[1] - Defina o numero de bits (ADC) para descrever o sinal")
    print("[2] - Defina a frequencia de amostragem [em Hz]")
    print("[3] - Rode a comparação")
    print("[4] - Sair")


print('\t\t\t Bem vindo! Este programa faz análise de quantização.')
print(
    '\t\t Este programa mostra um sinal seno continuo [0 a 5V], e os efeitos da amostragem e '
    'quantização escolhida no ADC.\n')

fc = 1
fd = 0
qt = 0

print_menu()
x = input()
option = int(x)

while option != 4:
    if option == 1:
        print("Digite o numero de bits do conversor Analógico-Digital")
        qt = int(input())
    elif option == 2:
        print("Digite a frequencia de amostragem [em Hz]")
        fd = int(input())
    elif option == 3:
        print("Rodando ...")
        plot_signals(freq_cont=fc, freq_disc=fd, quantization=qt)
    else:
        print("Escolha a opção correta!")
    print_menu()
    x = input()
    option = int(x)
