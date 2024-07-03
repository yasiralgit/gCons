import os
from subprocess import call
from Bio.Seq import Seq
from Bio import SeqIO
import sys
import csv
import time
import numpy
import matplotlib.pyplot as plt
import numpy as np
import mplcursors 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk
from customtkinter import * 


def affStructSeg(all_data_sets,tabFrame,argv,k,alpha,beta):
    num_colors = max(len(data_set) for data_set in all_data_sets) # on va avoir autant de couleur qu'il y a de groupes communs distincts
    cmap = plt.cm.get_cmap('viridis', num_colors) # on prend la bibliothèque de couleur viridis pour num_colors couleurs
    # pour l'échelle on récupère la plus petite et la plus grande valeur des zones communes
    #for elt in all_data_sets : print("\n\tELEMENT DE DATA SETS : ",elt)
    min_val = min(min(zone[0] for group in data_set for zone in group) for data_set in all_data_sets)
    max_val = max(max(zone[1] for group in data_set for zone in group) for data_set in all_data_sets)
    # on prend la taille du plus grand groupe commun de zone pour pouvoir espacé avec précaution chaque segment de chaque génome
    max_data_set = max( len(data_set) for data_set in all_data_sets)
    fig, ax = plt.subplots(figsize=(5,2))
    legend_handles = []
    #legend_handles2 = []
    segment_info = {}
    y_mid = []
    #tkinter n'accepte que les couleurs en hexa or on a des couleurs en codage RGB
    def rgba_to_hex(rgba):
        r, g, b, _ = rgba
        return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'

    labels = []
    for ds_index, data_set in enumerate(all_data_sets):
        y_value =3*max_data_set*ds_index # pour chaque génome on va avoir cet espacement
        y_mid.append(y_value+(len(data_set)/2)) # au milieu de cet espacement on va placer dans les ordonnées le nom du génome
        for group_index, group in enumerate(data_set):
            y_value += 1 # pour pouvoir écarter chaque zone pour éviter tout overlapping
            color = cmap(group_index/num_colors) #pour chaque zone d'un groupe on lui associe une couleur
            hex_color = rgba_to_hex(color) #on convertit cette couleur en hexa pour la phase tkinter ensuite
            if ds_index == 0 : # si on est dans le génome de référence
                # on récupère les bornes du groupe ce qui va nous servir de légénde pour le groupe 
                min_group = group[0][0] 
                max_group = group[-1][1]
                labels.append(f'{min_group}-{max_group}')
                #legend_handles2.append(plt.Line2D([0],[0],color=color,lw=4,label=f'{min_group}-{max_group}'))
                legend_handles.append((hex_color,f'{min_group}-{max_group}')) #pour la légende tkinter
            for zone in group:
                line, = ax.plot(zone, [y_value,y_value], linewidth=10, color=color) #on trace chaque zone
                #print(group_index)
                segment_info[line] = labels[group_index] #on associe à la zones les bornes du génome de référence
    ax.set_xlim(min_val - 1, max_val + 1) 
    ax.set_yticks(y_mid)
    ystickslables = []
    for fic_index,fic in enumerate(argv):
        ystickslables.append(f"Génome {fic_index+1} in {fic}") #on crée les modalités des ordonnées
    ax.set_yticklabels(ystickslables)
    fig.tight_layout(rect=[0, 0, 1, 1])

    def create_scrollable_legend(handles, title):
        root = tabFrame

        main_frame = ttk.Frame(root, padding='0.05i')

        name = "Reference genome = "+argv[0]+ "| compared genomes = "
        for fic in argv[1:] :
            name += fic + " "
        labelF = ttk.Label(main_frame, text= f"{name} ", background='white')
        labelF2 = ttk.Label(main_frame,text=f"| k-mer size = {k} | genome with the same k-mer = {alpha} | genome with the same link between k-mer = {beta}", background='white')

        #main_frame = CTk.Frame(root)
        labelF.pack(fill=tk.BOTH, expand=1)
        labelF2.pack(fill=tk.BOTH, expand=1)
 

        scroll_canvas = tk.Canvas(main_frame, background="blue")
        
        #Création de l'étiquette qui va gérer le clic sur zone
        info_label = tk.Label(main_frame, text="Cliquez sur un segment pour voir les informations du groupe", bg="white", fg="black", font=("Arial",12))
        info_label.pack(side=tk.BOTTOM, fill=tk.X) # cette étiquette est placée en bas de tout et prend toute la largeur

        scrollbarH = ttk.Scrollbar(main_frame, orient=tk.HORIZONTAL, command=scroll_canvas.xview)
        scrollbarH.pack(side=tk.BOTTOM, fill=tk.X)

        scrollbarV = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=scroll_canvas.yview)
        scrollbarV.pack(side=tk.RIGHT, fill=tk.Y)

        main_frame.pack(side=tk.TOP, expand=1, fill=tk.BOTH)
        scroll_canvas.pack(side=tk.LEFT,  expand=1, fill=tk.BOTH)

        scroll_canvas.configure(xscrollcommand=scrollbarH.set,yscrollcommand=scrollbarV.set)

        canvas_frame = ttk.Frame(scroll_canvas)
        #scroll_canvas.create_window((0,0), window=canvas_frame,anchor='center', height=300)
        canvas_frame.pack(expand=1,fill=tk.BOTH)
        canvas = FigureCanvasTkAgg(fig, master=canvas_frame) # on récupère la figure de matpolib pour l'afficher via tkinter
        canvas.draw()
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=1) # on place la figure matpotlib sur le coté gauche du widjet

        canvas_frame.update_idletasks()
        scroll_canvas.config(scrollregion=scroll_canvas.bbox("all"))

        #blocCan = tk.Canvas(canvas)
        #scrollbarCan = ttk.Scrollbar(main_frame,orient=tk.HORIZONTAL)
        #scrollbarCan.pack(side=tk.TOP, fill=tk.X)
#
        #scrollbarCan.config(command=canvas.xview)
        #main_frame.config(xscrollcommand=scrollbarCan.set)
        #main_frame.config(scrollregion=main_frame.bbox("all"))


        #ajout de la barre d'outil
        toolbar = NavigationToolbar2Tk(canvas, main_frame) #on cree une toolbar par rapport à canvas dans main_frame
        toolbar.update() #on met à jour la barre d'outil sur l'état actuel de la figure matpotlib
        #canvas_widget.pack( fill=tk.BOTH, expand=1) #on met le canvas au dessus de la toolbar dans main_frame


 
        legend_frame = ttk.Frame(main_frame) #on crée un widget de main_frame qui va etre la légende
        legend_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1) # 

        scrollbar = ttk.Scrollbar(legend_frame, orient=tk.VERTICAL) #on crée une scrollbar verticale du widget de la légende
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y) # on met la scrollbar à droite


        legend_canvas = tk.Canvas(legend_frame) #dans le widget da la légénde on crée un canvas qui va servir pour les icones de couleur
        legend_frame_inner = ttk.Frame(legend_canvas) #widget qui contient les éléments de la légende
        scrollbar.config(command=legend_canvas.yview) #on config la scrollbar pour le widget legend_canvas
        legend_canvas.create_window((0, 0), window=legend_frame_inner, anchor='nw') #crée une fenetre en haut à gauche de legend_canvas qui contient le widget legend_frame_inner
        legend_canvas.config(yscrollcommand=scrollbar.set) #on config de défilement de legend_canvas 
        legend_canvas.pack(side=tk.LEFT,fill=tk.BOTH, expand=1)

        for color,label in handles :
            label_frame = ttk.Frame(legend_frame_inner) #on crée un widget du widget d'élément de la légende
            color_label = tk.Label(label_frame, bg=color,width=2) # on crée un label de couleur du widget label_frame
            text_label = tk.Label(label_frame, text=label, anchor='w') # on crée un label de texte du widget label_frame ancrée à gauche
            color_label.pack(side=tk.LEFT, fill=tk.Y) # la couleur va prendre toute la hauteur disponible dans label_frame
            text_label.pack(side=tk.LEFT, fill=tk.BOTH, expand=1) # le texte va prendre toute la hauteur et la largeur disponible dans label_frame
            label_frame.pack(fill=tk.X, pady=1) #on place l'élement dans la légende on gardant un petit espace pour le suivant

        legend_frame_inner.update_idletasks() #mets à jours dans tkinter les modifications qu'on vient d'effectuer
        legend_canvas.config(scrollregion=legend_canvas.bbox("all")) # configure tout le widget legend_canvas comme étant défilable 


        def on_click(event): 
            x, y = event.xdata, event.ydata
            if x is not None and y is not None : #on vérifie que les coordonées du clic sont correctes
                for line in segment_info : #pour chaque zone on va voir si on en a pas une qui coordonne avec ce clic
                    line_xdata = line.get_xdata() # on récupère les abscisses 
                    line_ydata = line.get_ydata() # on récupère les ordonnées
                    if line_xdata[0] <= x <= line_xdata[1] and line_ydata[0] - 1 <= y <= line_ydata[0] + 1 : #verifie que le clic est dans la zone du segment 
                        group_id = segment_info[line] # donne le label de cette zone/ de ce segment 
                        print(f"Clicked on {group_id}")
                        info_label.config(text=f"Vous avez cliqué sur la zone : {group_id}")
                        break # si une ligne correspond on a plus besoin de continuer d'itérer
        canvas.mpl_connect("button_press_event",on_click) #appelle on_click lors d'un clic de souris dans le widget qui contient la figure matpotlib        
        #canvas.mpl_connect("scroll_event",on_scroll)


        scroll_canvas.update_idletasks()
        scroll_canvas.config(scrollregion=scroll_canvas.bbox("all"))

    create_scrollable_legend(legend_handles,"Group IDs")

#Pour debogger
def affichMatGen(matGen):
    print("[")
    for i in range (len(matGen)):
        print(f"\t genome {i+1}")
        print(f"\t [")
        for j in range (len(matGen[i])) :
            print(f"\t {matGen[i][j]}")
        print(f"\t ]")
    print("]")

def affichMatGenComp(matGen,ind):
    print("[")
    for i in range (len(matGen)):
        print(f"\t genome {i+1} zone {ind}")
        print(f"\t [")
        print(f"\t {matGen[i][ind]}")
        print(f"\t taille de la zone => {len(matGen[i][ind])}")
        print(f"\t ]")
    print("]")

def affichMatGenCompMeilleur(matGen) :
    for i in range (len(matGen[0])):
        affichMatGenComp(matGen,i)

def affichMatGenRef(matGen):
    print("[")
    print("\t genome 1")
    print("\t [")
    for j in range (len(matGen[0])) :
        print(f"\t {matGen[0][j]}")
    print("\t ]")
    print("]")

def appelKmer():
    k = int(input("Saisissez la taille de kmer souhaitée : "))
    return k
''' Finalement gamma pas utilisé
def calculGamma(g,Tposition,saveGkampi):
    nb = 0 
    for i in range(Tposition[g]) :
        nb += len(Tposition[g][i])
    return nb/len(saveGkampi[g])
'''






def lancerGkampi(k,fichier,num):
    """
    Permeet de lancer Gkampi pour un fichier passé en entré
    """
    commandeGkampi = ["./gkampi"]
    commandeGkampi.append(fichier)
    commandeGkampi.append("--output")
    commandeGkampi.append("result"+str(num)+".csv")
    commandeGkampi.append("--pos")
    commandeGkampi.append("--column")
    commandeGkampi.append("--kmer-len")
    commandeGkampi.append(str(k))
    commandeGkampi.append("--quiet")
    call(commandeGkampi)
    print (commandeGkampi)
    return 0


def lancerRedOak(k,argv):
    commandeRedOak = ["./redoak"]
    for element in argv:
        commandeRedOak.append("--genome")
        commandeRedOak.append(element)
    commandeRedOak.append("--output")
    commandeRedOak.append("test.txt")
    commandeRedOak.append("--kmer")
    commandeRedOak.append(str(k))
    print(commandeRedOak)
    call(commandeRedOak)
    return 0 


def commandeGrep(alpha,argv):
    commandeGrep = "grep -E \""
    val = int(alpha)
    while val < len(argv)+1:
        commandeGrep = commandeGrep + "\("+str(val)+"\)"
        val = val + 1
        if val < len(argv)+1:
            commandeGrep=commandeGrep+"|"
    commandeGrep = commandeGrep + "\" test.txt > commun.txt" 
    print(commandeGrep)
    os.system(commandeGrep)


def listeKmers():
    liste = []
    with open("commun.txt","r") as File:
        cpt = 0
        for ligne in File:
            lettre = 0
            k_mers = ""
            while ligne[lettre] != " " and (ligne[lettre] == "A" or ligne[lettre] == "C" or ligne[lettre] == "T" or ligne[lettre] == "G"):
                k_mers = k_mers + ligne[lettre]
                lettre = lettre+1
            liste.append(k_mers) 
    print("Nombre de k-mers identiques == > "+ str(len(liste)))
    return liste 

def listePosition(k,argv):
    liste = listeKmers()
    Tposition=[]
    saveGkampi = {}
    saveGkampi2 = {}
    num = 1
    for element in argv:
        print(argv)
        cpt = 1
        lancerGkampi(k,element,num)
        saveGkampi["g"+str(num)]= {}
        saveGkampi2["g"+str(num)]= {}
        Tposition.append([]) 
        for i in range(len(liste)):
            Tposition[-1].append([])
        debut = time.time()
        with open("result"+str(num)+".csv",newline='') as File: 
            Lecture = csv.reader(File,delimiter=" ") 
            Lecture = list(Lecture) 
            if cpt == 1:
                for i in Lecture:
                    saveGkampi["g"+str(num)][int(i[1][1:])] = i[0] 
                    saveGkampi2["g"+str(num)][i[0]] = int(i[1][1:]) #partiel
            print("La taille est de : ",len(liste))
            for i in range (len(liste)):
                gkampIndex = 0
                while gkampIndex<len(Lecture) and str(Lecture[gkampIndex][0]) != str(liste[i])  : 
                    gkampIndex = gkampIndex +1 
                while gkampIndex<len(Lecture) and str(Lecture[gkampIndex][0]) == str(liste[i]) :
                    Tposition[-1][i].append(int(Lecture[gkampIndex][1][1:]))
                    gkampIndex = gkampIndex + 1          
        cpt += 1
        num += 1 
        print(time.time()-debut)
    return Tposition,saveGkampi,saveGkampi2


def extensible(pos,NKmer,Tposition,beta): 
    """
    Permet de déterminer si un K-mers peu être étendue et renvoie si il peut être étendue et aussi comment ce k-mer peut être étendue
    """
    tabRes = []
    nbtrouve = 0
    for Kmer in range (0,len(Tposition[0])):
            for position in range (0,len(Tposition[0][Kmer])):
                if Tposition[0][Kmer][position] == (pos+1):
                    ksauv = Kmer
                    tabRes.append([Tposition[0][Kmer][position]])
                    nbtrouve = 1
                    genomeVu = []
                    for genome in range(1,len(Tposition)):
                        tabRes.append([])
                        for position2 in range (0,len(Tposition[genome][NKmer])):
                            if (Tposition[genome][NKmer][position2]+1) in Tposition[genome][Kmer]:
                                if genome not in genomeVu :
                                    tabRes[-1].append(Tposition[genome][NKmer][position2]+1)
                                    nbtrouve += 1
                                    genomeVu.append(genome)

    #print("\n\tNB TROUVE : ",nbtrouve,"et beta : ",beta)
    if (nbtrouve>=beta and nbtrouve != 0):
        return True,ksauv,tabRes
    else :
        return False,-1,-1

def simpCent(matGenZone):
    matGenZoneFinale = [] #initialement vide
    for j in range (len(matGenZone)):
        matGenZoneFinale.append([]) #aucun groupe pour l'instant
    for groupe_index,groupe in enumerate(matGenZone[0]): 
        if groupe != [] :
            bool = True
            for groupeP in range(len(matGenZoneFinale[0])): 
                for zone in groupe: 
                    if zone in matGenZoneFinale[0][groupeP]: 
                        if len(groupe) > len(matGenZoneFinale[0][groupeP]): 
                            matGenZoneFinale[0][groupeP] = groupe 
                            bool = False
                            for i in range(1,len(matGenZone)):
                                matGenZoneFinale[i][groupeP] = matGenZone[i][groupe_index]
            if bool == True:
                matGenZoneFinale[0].append(groupe)
                for i in range (1,len(matGenZone)):
                    matGenZoneFinale[i].append(groupe)
    for elt in matGenZoneFinale :
        elt = elt.sort
    return matGenZoneFinale

def triCent(matGenZone):
    matGenRes = []
    for i in range (len(matGenZone)) :
        matGenRes.append([])
    clas = []
    recup = {}
    for groupe_index,groupe in enumerate(matGenZone[0]) :
        recup[groupe[0][0]] = groupe_index
        clas.append(groupe[0][0])
    #print(clas)
    clas.sort()
    for min in clas :
        for i in range(len(matGenZone)) :
            matGenRes[i].append(matGenZone[i][recup[min]])
    #print(matGenRes)
    return matGenRes
            

def zonesCommunes(Tposition,beta):
    liste_Finale = []
    listeDesParcourus = []
    tabCentAdr = []
    for i in range (len(Tposition)) :
        tabCentAdr.append([]) #un tab pour chaque génome
    for i in range (len(Tposition[0])): #pour chaque kmer du génome de référence
        res = []
        for z in range(len(Tposition[0][i])):#pour chaque position de ce kmer on va étendre au maximum si cette position n'a pas déjà été vu
            cpt = 1 
            if Tposition[0][i][z] not in listeDesParcourus :
                listeDesParcourus.append(Tposition[0][i][z])
                res.append([Tposition[0][i][z]])
                for q in range (len(Tposition)): #pour chaque kmer commun du génome de référence on commence à créer une zone dans les autres génomes, elle restera vide si il n y a pas de liaison qui satisfasse beta 
                    tabCentAdr[q].append([])
                bool,kmer,tabRes = extensible(Tposition[0][i][z],i,Tposition,beta)
                #print(f"\n\{Tposition[0][i][z]} : ",Tposition[0][i][z])
                #print("\n\tBOOL : ",bool)
                while bool:
                    for q in range (len(tabRes)):
                        if tabRes[q] != [] :
                            tabCentAdr[q][-1].append([tabRes[q][0]-1,tabRes[q][0]]) #on ajoute à la zone qui vient d'etre crée la liaison
                    res[-1].append(tabRes[0][0])
                    listeDesParcourus.append(tabRes[0][0])
                    bool,kmer,tabRes = extensible(tabRes[0][0],kmer,Tposition,beta)
        #Pour ne garder que les plus grandes zones distinctes
        for element in res: 
            bool = True
            for element2 in range(len(liste_Finale)): 
                for pos in element: 
                    if (pos) in liste_Finale[element2]: 
                        if len(element) > len(liste_Finale[element2]): 
                            liste_Finale[element2] = element 
                            bool = False
            if bool == True:
                liste_Finale.append(element)
    liste_Finale.sort()
    #print("liste_Finale : ",liste_Finale)
    return liste_Finale,tabCentAdr#,tabCentAdrSimp

def adr2seq(liste_Finale,saveGkampi):
    listeseq = []
    for i in liste_Finale:
        sequence = saveGkampi["g1"][i[0]] 
        for y in i[1:]:
            sequence = sequence + saveGkampi["g1"][y][-1] 
        listeseq.append(sequence)
    return listeseq


def sauvGenCons(liste_Finale,listeseq,k,argv):
    scaffoldList = []
    nameliste = []
    infoSup=""
    with open("ResulFasta.fasta","w") as File:
        lecture = SeqIO.to_dict(SeqIO.parse(argv[0],"fasta")) 
        for seq in lecture:
            if scaffoldList == []:
                scaffoldList.append([0,len(lecture[seq])])
                nameliste.append(seq)
            else :
                scaffoldList.append([(scaffoldList[-1][1])+1,len(lecture[seq])+scaffoldList[-1][1]])
                nameliste.append(seq)
        for id,sequence in enumerate(listeseq):
            for compteur,pos in enumerate(scaffoldList):
                infoSup = f" | {liste_Finale[id][0]},{liste_Finale[id][-1] + k-1} VS {pos[0]},{pos[1]} "
                if liste_Finale[id][0]>= pos[0] and liste_Finale[id][-1] + k-1 <= pos[1]: # si il est contenu dans le scaffold
                    infoSup += "| Pos in seq ("+ nameliste[compteur] +") inside the reference fasta file : begin = " + str( -pos[0] + liste_Finale[id][0]) 
            entete = ">sequence/scaffold_"+ str(id+1) +" | Position with Gkampi index in g1 "+ str(liste_Finale[id][0]) +" ==> "+ str(liste_Finale[id][-1] + k-1) + " | " + str(argv[0]) + infoSup +"\n"
            File.write(entete+sequence+"\n")
    return 0 

def sauvZonesCommunes(liste_Finale):
    with open("resultat.csv",'w',newline='') as File: # On écrit nos liste de positions dans un fichier CSV pour pouvoir faire des vérifications
        Ecriture = csv.writer(File,delimiter=";")
        Ecriture.writerows(liste_Finale)


def main(argv,kmer,alpha,beta,tabFrame) :
    alpha = int(numpy.round((int(alpha) * len(argv)) / 100)) 
    beta = int(numpy.round((int(beta) * len(argv)) / 100)) 
    kmer = int(kmer)
    lancerRedOak(kmer,argv)
    commandeGrep(alpha,argv)
    listeKmersCommuns = listeKmers()
    listeKmersPositions,dicoAdr2Kmer,dicoKmer2Adr = listePosition(kmer,argv)
    #for elt in listeKmersPositions : print("\n\tELEMENT DE TPOSITION : ",elt)
    listeZonesCommunes,data_sets = zonesCommunes(listeKmersPositions,beta)
    #for elt in data_sets : print("\n\tELEMENT DE DATASET : ",elt)
    affichMatGenCompMeilleur(triCent(simpCent(data_sets)))

    listeSeq = adr2seq(listeZonesCommunes,dicoAdr2Kmer)
    #print(listeKmersPositions)
    sauvGenCons(listeZonesCommunes,listeSeq,kmer,argv)
    sauvZonesCommunes(listeZonesCommunes)


    affStructSeg(triCent(simpCent(data_sets)),tabFrame,argv,kmer,alpha,beta)