import os
from subprocess import call
from Bio import SeqIO
import csv
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
from tkinter import ttk
from customtkinter import * 

nbSub = 0 #nombre de Submit envoyé

#-------------Pour debugger----------------------------------------------------------------------------------
def affichMatGen(matGen):
    '''
    Affiche dans une structure de donnée adaptée tous les éléments contenus dans chaque génome
    '''
    print("[")
    for i in range (len(matGen)):
        print(f"\t genome {i+1}")
        print(f"\t [")
        for j in range (len(matGen[i])) :
            print(f"\t {matGen[i][j]}")
        print(f"\t ]")
    print("]")

def affichMatGenComp(matGen,ind):
    '''
    Affiche dans une structure de donnée adaptée l'élément d'indice indiqué contenu dans chaque génome
    '''
    print("[")
    for i in range (len(matGen)):
        print(f"\t genome {i+1} zone {ind}")
        print(f"\t [")
        print(f"\t {matGen[i][ind]}")
        print(f"\t taille de la zone => {len(matGen[i][ind])}")
        print(f"\t ]")
    print("]")

def affichMatGenCompMeilleur(matGen) :
    '''
    Affiche dans une structure de donnée adaptée chaque élément de chaque génome en simultané de manière à pouvoir les comparer
    '''
    for i in range (len(matGen[0])):
        affichMatGenComp(matGen,i)

def affichMatGenRef(matGen):
    '''
    Affiche dans une structure de donnée adaptée tous les éléments du génome de référence
    '''
    print("[")
    print("\t genome 1")
    print("\t [")
    for j in range (len(matGen[0])) :
        print(f"\t {matGen[0][j]}")
    print("\t ]")
    print("]")
#-------------fin -------------------------------------------------------------------------------------------


def lancerGkampi(k,fichier,num):
    '''
    Permet de lancer gkampi pour le fichier et la taille de k-mer indiqués, num permet de nommer le fichier de sortie
    \param k : entier
    \param fichier : string
    \param num : entier
    \return 0 
    '''
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
    '''
    Permet de lancer redOak pour la taille de k-mer les fichiers indiqués
    \param k : entier
    \param argv : tableau de string
    \return 0
    '''
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
    '''
    Permet de lancer la commande grep qui va filtrer les k-mers communs du fichier de sortie redOak selon alpha
    \param alpha : entier
    \param argv : tableau de string
    '''
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


def fliste():
    '''
    Crée la structure de donnée "liste" qui contient tous les k-mers communs à partir du fichier les contenant "commun.txt"
    \return liste : tableau de string
    '''
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

def fTposition(k,argv):
    '''
    Crée la structure de donnée Tposition qui contient toutes les positions des k-mers communs dans chaque génome
    \param k : entier
    \param argv : tableau de string
    \return Tposition : tableau de tableau de tableau d'entiers
    \return saveGkampi : dictionnaire qui associe à chaque position du génome de référence le k-mer associé
    '''
    liste = fliste()
    Tposition=[]
    saveGkampi = {}
    num = 1
    for element in argv:
        print(argv)
        lancerGkampi(k,element,num)
        Tposition.append([])
        for i in range(len(liste)):
            Tposition[-1].append([])
        debut = time.time()
        with open("result"+str(num)+".csv",newline='') as File: 
            Lecture = csv.reader(File,delimiter=" ") 
            Lecture = list(Lecture) 
            if num == 1:
                for i in Lecture:
                    saveGkampi[int(i[1][1:])] = i[0] 
            print("La taille est de : ",len(liste))
            for i in range (len(liste)):
                gkampIndex = 0
                while gkampIndex<len(Lecture) and str(Lecture[gkampIndex][0]) != str(liste[i])  : 
                    gkampIndex = gkampIndex +1 
                while gkampIndex<len(Lecture) and str(Lecture[gkampIndex][0]) == str(liste[i]) :
                    Tposition[-1][i].append(int(Lecture[gkampIndex][1][1:]))
                    gkampIndex = gkampIndex + 1          
        num += 1 
        print(time.time()-debut)
    return Tposition,saveGkampi


def extensible(pos,NKmer,Tposition,beta): 
    '''
    Permet de déterminer si un k-mer peut être étendu
    \param pos : entier
    \param NKmer : entier
    \param Tposition : tableau de tableau de tableau d'entiers
    \param beta : entier
    \return True l'indice du k-mer associé et sa position si le k-mer peut etre étendu False,-1,-1 sinon
    '''
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

    if (nbtrouve>=beta and nbtrouve != 0):
        return True,ksauv,tabRes
    else :
        return False,-1,-1

def simpCent(matGenZone):
    '''
    Réduit la structure de donnée matGenZone en enlevant les doublons
    \param matGenZone : tableau de tableau de tableau de tableau d'entiers
    \return matGenZoneFinale : tableau de tableau de tableau de tableau d'entiers
    '''
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
    '''
    Trie la structure de donnée matGenZone par ordre croissant selon le génome de référence
    \param matGenZone : tableau de tableau de tableau de tableau d'entiers
    \return matGenRes : tableau de tableau de tableau de tableau d'entiers
    '''
    matGenRes = []
    for i in range (len(matGenZone)) :
        matGenRes.append([])
    clas = []
    recup = {}
    for groupe_index,groupe in enumerate(matGenZone[0]) :
        recup[groupe[0][0]] = groupe_index
        clas.append(groupe[0][0])
    clas.sort()
    for min in clas :
        for i in range(len(matGenZone)) :
            matGenRes[i].append(matGenZone[i][recup[min]])
    return matGenRes
            

def zonesCommunes(Tposition,beta):
    '''
    Crée la structure de donnée liste_Finale qui contient la succession de sommets constituant les zones communes dans le génome de référence ainsi que tabCentAdr qui comporte les zones communes pour chaque génome en sucession d'arête
    \param Tposition : tableau de tableau de tableau d'entiers
    \param beta : entier
    \return liste_Finale : tableau d'entiers
    \return tabCentAdr : tableau de tableau de tableau de tableau d'entiers
    '''
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
                while bool:
                    for q in range (len(tabRes)):
                        if tabRes[q] != [] :
                            tabCentAdr[q][-1].append([tabRes[q][0]-1,tabRes[q][0]]) #on ajoute à la zone qui vient d'etre crée la liaison
                    res[-1].append(tabRes[0][0])
                    listeDesParcourus.append(tabRes[0][0])
                    bool,kmer,tabRes = extensible(tabRes[0][0],kmer,Tposition,beta)
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
    return liste_Finale,tabCentAdr#,tabCentAdrSimp

def adr2seq(liste_Finale,saveGkampi):
    '''
    Convertit une liste de succession de positions du génome de référence en liste de séquence de nucélotides
    \param liste_Finale : tableau d'entiers
    \param saveGkampi : dictionnaire
    \return listeseq : tableau de string
    '''
    listeseq = []
    for i in liste_Finale:
        sequence = saveGkampi[i[0]] 
        for y in i[1:]:
            sequence = sequence + saveGkampi[y][-1] 
        listeseq.append(sequence)
    return listeseq


def sauvGenCons(liste_Finale,listeseq,k,argv):
    '''
    Crée le fichier de sortie fasta qui contient le génome consensus
    \param liste_Finale : tableau d'entiers
    \param listeseq : tableau de string
    \param k : entier
    \param argv : tableau de string
    \return 0 
    '''
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
    '''
    Crée le fichier de sortie resultat.csv qui contient les zones communes pour permettre le débuggage
    \param liste_Finale : tableau d'entiers
    \return 0
    '''
    with open("resultat.csv",'w',newline='') as File: # On écrit nos liste de positions dans un fichier CSV pour pouvoir faire des vérifications
        Ecriture = csv.writer(File,delimiter=";")
        Ecriture.writerows(liste_Finale)
    return 0 




def affStructSeg(tabCentAdr,tabFrame,app,argv,k,alpha,beta,nbSub):
    '''
    Crée et affiche la figure matplotlib des zones communes
    \param tabCentAdr : tableau de tableau de tableau de tableau d'entiers, matrice des zones communes pour chaque génome
    \param tabFrame : objet customtkinter CTkScrollableFrame
    \param app : objet racine de l'interface customtkinter
    \param argv : tableau de string, noms des fichiers insérés
    \param k : entier, taille de k-mer inséré
    \param alpha : entier, nombre de génome validant alpha
    \param beta : entier, nombre de génome validant béta
    \param nbSub : entier, nombre de submit envoyé
    \return 0 
    '''
    #configuration des couleurs associés aux zones communes via la bibliothèque viridis
    num_colors = max(len(data_set) for data_set in tabCentAdr) 
    cmap = plt.cm.get_cmap('viridis', num_colors)
    # pour l'échelle du schéma on récupère la plus petite et la plus grande position 
    min_val = min(min(zone[0] for group in data_set for zone in group) for data_set in tabCentAdr)
    max_val = max(max(zone[1] for group in data_set for zone in group) for data_set in tabCentAdr)
    # on prend la taille du plus grand groupe de zone pour pouvoir espacé avec précaution en ordonnée chaque segment de chaque génome
    max_data_set = max( len(data_set) for data_set in tabCentAdr)
    #on crée le schéma
    fig, ax = plt.subplots(figsize=(5,2))
    legend_handles = []
    segment_info = {}
    y_mid = []
    #tkinter n'accepte que les couleurs en codage hexadécimal or on a des couleurs en codage RGB, une fonction de conversion est donc nécessaire
    def rgba_to_hex(rgba):
        '''
        \param rgba codage de couleur en rgb
        \return codage de couleur en hexadécimal
        '''
        r, g, b, _ = rgba
        return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'

    labels = []
    for ds_index, data_set in enumerate(tabCentAdr):
        y_value =3*max_data_set*ds_index # pour chaque génome on va avoir cet espacement
        y_mid.append(y_value+(len(data_set)/2)) # au milieu de cet espacement on va placer dans les ordonnées le génome suivi de son numéro
        for group_index, group in enumerate(data_set):
            y_value += 1 # pour pouvoir écarter chaque groupe de zones pour éviter tout chevauchement
            color = cmap(group_index/num_colors) #pour chaque groupe on associe une couleur
            hex_color = rgba_to_hex(color) #on convertit cette couleur en hexadécimal
            if ds_index == 0 : 
                # on récupère les bornes du groupe dans le génome de référence ce qui va nous servir de légénde
                min_group = group[0][0] 
                max_group = group[-1][1]
                labels.append(f'{min_group}-{max_group}')
                legend_handles.append((hex_color,f'{min_group}-{max_group+k-1}')) 
            for zone in group:
                #on trace chaque zone
                line, = ax.plot(zone, [y_value,y_value], linewidth=10, color=color) 
                segment_info[line] = labels[group_index] #on associe à chaque zone la légende de son groupe
    ax.set_xlim(min_val - 1, max_val + 1) 
    ax.set_yticks(y_mid)
    ystickslables = []
    for fic_index,fic in enumerate(argv):
        ystickslables.append(f"Génome {fic_index+1} in {fic}") 
    ax.set_yticklabels(ystickslables)
    fig.tight_layout(rect=[0, 0, 1, 1])
    def on_closing():
        '''
        Permet de bien fermer la figure matplotlib quand on ferme l'interface customtkinter
        '''
        plt.close('all')  # Fermer toutes les figures matplotlib
        after_tasks = app.tk.eval("after info").split()
        for after_id in after_tasks: # Fermer toutes les fonctionnalités associées à l'objet app
            app.after_cancel(after_id)
        print("fermeture")
        app.destroy()  # Détruire l'application tkinter
    app.protocol("WM_DELETE_WINDOW", on_closing)

    def create_scrollable_legend(handles,indice):
        '''
        Met en place la fenêtre tkinter qui va heberger la figure matplotlib ainsi que sa légende
        '''
        print(f"Nombre de schéma : {indice}")
        root = tabFrame
        main_frame = ttk.Frame(root, padding='0.05i', name=f"mf{indice}")

        #Création des labels d'information qui renseigne les paramètres du schéma affiché
        name = "Reference genome = "+argv[0]+ "| compared genomes = "
        for fic in argv[1:] :
            name += fic + " "
        labelF = ttk.Label(main_frame, text= f"{name} ", background='white')
        labelF2 = ttk.Label(main_frame,text=f"| k-mer size = {k} | genome with the same k-mer = {alpha} | genome with the same link between k-mer = {beta}", background='white')
        labelF.pack(fill=tk.BOTH, expand=1)
        labelF2.pack(fill=tk.BOTH, expand=1)
 
        
        #Création de l'étiquette qui va gérer le clic sur zone
        info_label = tk.Label(main_frame, text="Cliquez sur un segment pour voir les informations du groupe", bg="white", fg="black", font=("Arial",12))
        info_label.pack(side=tk.BOTTOM, fill=tk.X) # cette étiquette est placée en bas de tout et prend toute la largeur

        #Création d'un objet Canva qui va contenir la figure matplotlib
        scroll_canvas = tk.Canvas(main_frame, background="blue")
        scrollbarH = ttk.Scrollbar(main_frame, orient=tk.HORIZONTAL, command=scroll_canvas.xview)
        scrollbarH.pack(side=tk.BOTTOM, fill=tk.X)
        scrollbarV = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=scroll_canvas.yview)
        scrollbarV.pack(side=tk.RIGHT, fill=tk.Y)
        main_frame.pack(side=tk.TOP, expand=1, fill=tk.BOTH)
        scroll_canvas.pack(side=tk.LEFT,  expand=1, fill=tk.BOTH)
        scroll_canvas.configure(xscrollcommand=scrollbarH.set,yscrollcommand=scrollbarV.set)
        canvas_frame = ttk.Frame(scroll_canvas)
        scroll_canvas.create_window((0,0), window=canvas_frame,anchor='center', height=300)
        canvas_frame.pack(expand=1,fill=tk.BOTH)
        canvas = FigureCanvasTkAgg(fig, master=canvas_frame) # on récupère la figure de matpolib pour l'afficher via tkinter
        canvas.draw()
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=1) # on place la figure matpotlib sur le coté gauche du widjet
        canvas_frame.update_idletasks()
        scroll_canvas.config(scrollregion=scroll_canvas.bbox("all"))

        #Ajout de la Toolbar de matplotlib en widget tkinter
        toolbar = NavigationToolbar2Tk(canvas, main_frame) #on cree une toolbar par rapport à canvas dans main_frame
        toolbar.update() #on met à jour la barre d'outil sur l'état actuel de la figure matpotlib
        canvas_widget.pack( fill=tk.BOTH, expand=1) #on met le canvas au dessus de la toolbar dans main_frame


        #Création de l'objet Frame de la légende du schéma
        legend_frame = ttk.Frame(main_frame) 
        legend_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1) 
        scrollbar = ttk.Scrollbar(legend_frame, orient=tk.VERTICAL) #on crée une scrollbar verticale pour legend_frame
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y) # on place la scrollbar à droite
        legend_canvas = tk.Canvas(legend_frame) #création d'un objet Canvas à partir de legend_frame
        legend_frame_inner = ttk.Frame(legend_canvas) #création d'un widget pour chaque élément de la légende
        scrollbar.config(command=legend_canvas.yview) 
        legend_canvas.create_window((0, 0), window=legend_frame_inner, anchor='nw') #placement de legend_frame dans une fenêtre en haut à gauche du Canvas
        legend_canvas.config(yscrollcommand=scrollbar.set)
        legend_canvas.pack(side=tk.LEFT,fill=tk.BOTH, expand=1)

        #Création de chaque élément de la légende à partir de la description et de la couleur de chaque groupe de zone
        for color,label in handles :
            label_frame = ttk.Frame(legend_frame_inner) 
            color_label = tk.Label(label_frame, bg=color,width=2) 
            text_label = tk.Label(label_frame, text=label, anchor='w') 
            color_label.pack(side=tk.LEFT, fill=tk.Y) 
            text_label.pack(side=tk.LEFT, fill=tk.BOTH, expand=1) 
            label_frame.pack(fill=tk.X, pady=1) 

        legend_frame_inner.update_idletasks() #mise à jours dans tkinter des modifications effectuées
        legend_canvas.config(scrollregion=legend_canvas.bbox("all"))

        def on_click(event): 
            '''
            configure l'affichage de la zone d'un segment lorsqu'il est cliqué
            '''
            x, y = event.xdata, event.ydata
            if x is not None and y is not None : 
                for line in segment_info : 
                    line_xdata = line.get_xdata() 
                    line_ydata = line.get_ydata() 
                    if line_xdata[0] <= x <= line_xdata[1] and line_ydata[0] - 1 <= y <= line_ydata[0] + 1 : #la zone de clic en ordonnée est de 2
                        group_id = segment_info[line] 
                        print(f"Clicked on {group_id}")
                        info_label.config(text=f"Vous avez cliqué sur la zone : {group_id}")
                        break 
        canvas.mpl_connect("button_press_event",on_click)        

        def sup() :
            '''
            configure la suppression de schéma via le bouton associé en fonction du nombre de schéma
            '''
            global nbSub
            root.nametowidget(f"mf{nbSub}").destroy() #Destruction de l'objet Frame qui contient le schéma
            if nbSub == 1 :
                root.nametowidget(f"btn1").destroy() #Si il ne reste que le premier schéma : destruction du bouton de suppression
            nbSub -= 1             
        scroll_canvas.update_idletasks()
        scroll_canvas.config(scrollregion=scroll_canvas.bbox("all"))
        if indice == 1 :
            btn1 = ttk.Button(root, text ="Supprimer le dernier schéma", command = sup, name= f"btn{indice}")#Création du bouton de suppression de schéma au premier submit
            btn1.pack(pady = 10) 
    create_scrollable_legend(legend_handles,nbSub)
    return 0
    


def main(argv,kmer,alpha,beta,tabFrame2,app) :
    '''
    La fonction qui lance le programme
    \param argv : tableau de string 
    \param kmer : string
    \param alpha : string
    \param beta : string
    \param tabFrame2 : objet customtkinter CTkScrollableFrame
    \param app : objet racine de l'interface customtkinter
    '''
    global nbSub  
    nbSub+=1

    #Traitement des entrées et appel de redOak et gkampi
    alpha = int(numpy.round((int(alpha) * len(argv)) / 100)) 
    beta = int(numpy.round((int(beta) * len(argv)) / 100))
    kmer = int(kmer)
    lancerRedOak(kmer,argv)
    commandeGrep(alpha,argv)
    Tposition,saveGKampi = fTposition(kmer,argv)

    #Création des zones communes
    liste_Finale,tabCentAdr = zonesCommunes(Tposition,beta)
    affichMatGenCompMeilleur(triCent(simpCent(tabCentAdr)))
    
    #for key,value in saveGKampi.items() : print(f"\n\t{key}->{value}")
    listeSeq = adr2seq(liste_Finale,saveGKampi)

    #Les sorties finales
    sauvGenCons(liste_Finale,listeSeq,kmer,argv) 
    sauvZonesCommunes(liste_Finale)

    #Affichage du schéma des zones communes
    affStructSeg(triCent(simpCent(tabCentAdr)),tabFrame2,app,argv,kmer,alpha,beta,nbSub)

