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

nbSub = 0 #variable globale, nombre de Submit envoyé

#-------------Pour déboger----------------------------------------------------------------------------------
def affichMatGen(matGen):
    '''
    Affiche dans une structure de données adaptée tous les éléments contenus dans chaque génome.
    \param matGen : tableau de tableaux d'éléments
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
    Affiche dans une structure de données adaptée l'élément d'indice indiqué contenu dans chaque génome.
    \param matGen : tableau de tableaux d'éléments
    \param ind : entier
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
    Affiche dans une structure de données adaptée chaque élément de chaque génome simultanément, de manière à pouvoir les comparer.
    \param matGen : tableau de tableaux d'éléments
    '''
    for i in range (len(matGen[0])):
        affichMatGenComp(matGen,i)

def affichMatGenRef(matGen):
    '''
    Affiche dans une structure de données adaptée tous les éléments du génome de référence.
    \param matGen : tableau de tableaux d'éléments
    '''
    print("[")
    print("\t genome 1")
    print("\t [")
    for j in range (len(matGen[0])) :
        print(f"\n\t zone {j}")
        print(f"\t [")
        print(f"\t {matGen[0][j]}")
        print(f"\t ]")
        print(f"\t taille de la zone => {len(matGen[0][j])}")
    print("\t ]")
    print("]")
#-------------fin -------------------------------------------------------------------------------------------


def lancerGkampi(k,fichier,num):
    '''
    Permet de lancer gkampi pour le fichier et la taille de k-mer indiqués. \param num permet de nommer le fichier de sortie.
    \param k : entier, taille de k-mer
    \param fichier : chaîne de caractères, nom du fichier
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
    Permet de lancer redOak pour la taille de k-mer et les fichiers indiqués.
    \param k : entier
    \param argv : tableau de chaînes de caractères, tableau contenant le nom des fichiers
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
    Permet de lancer la commande grep qui va filtrer les k-mers communs du fichier de sortie redOak selon \param alpha.
    \param alpha : entier, seuil de génome à partir duquel on considère un k-mer comme commun
    \param argv : tableau de chaînes de caractères
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
    Crée la structure de données "liste" qui contient tous les k-mers communs à partir du fichier "commun.txt".
    \return liste : tableau de chaînes de caractères, tableau des k-mers communs
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
    Crée la structure de données \return Tposition qui contient toutes les positions des k-mers communs dans chaque génome.
    \param k : entier
    \param argv : tableau de chaînes de caractères
    \return Tposition : tableau de tableaux de tableaux d'entiers, tableau des positions des k-mers communs pour chaque génome
    \return saveGkampi : dictionnaire qui associe à chaque position du génome de référence le k-mer correspondant
    '''
    liste = fliste()
    Tposition=[]
    saveGkampi = {}
    num = 1
    for fichier in argv: #pour chaque génome
        lancerGkampi(k,fichier,num)
        Tposition.append([])
        #Ajout d'autant de tableau qu'il y a de k-mers communs dans liste
        for i in range(len(liste)):
            Tposition[-1].append([])
        debut = time.time()
        #Ouverture du fichier de sortie gkampi
        with open("result"+str(num)+".csv",newline='') as File: 
            Lecture = csv.reader(File,delimiter=" ") 
            Lecture = list(Lecture) 
            if num == 1:
                #S'il s'agit du génome de référence on alimente le dictionnaire saveGkampi
                for i in Lecture:
                    saveGkampi[int(i[1][1:])] = i[0] 
            print("La taille est de : ",len(liste))
            gkampIndex = 0 #Indice qui va nous permettre de d'itérer sur les lignes du fichier (*)
            for i in range (len(liste)): #pour chaque k-mer commun
                while gkampIndex<len(Lecture) and str(Lecture[gkampIndex][0]) != str(liste[i])  :#tant que le k-mer de la ligne n'est pas équivalent au k-mer commun
                    gkampIndex = gkampIndex +1 
                while gkampIndex<len(Lecture) and str(Lecture[gkampIndex][0]) == str(liste[i]) :#tant que le k-mer de la ligne est équivalent au k-mer commun
                    Tposition[-1][i].append(int(Lecture[gkampIndex][1][1:])) #Ajout de sa position dans Tposition
                    gkampIndex = gkampIndex + 1          
        num += 1 
        print(time.time()-debut)
    return Tposition,saveGkampi
# (*) On initialise gkampIndex avant la boucle de liste car redOak et gkampi indexent de même manière les génomes, il n y a pas besoin de remettre gkampIndex à 0 pour chaque k-mer commun.


def extensible(pos,NKmer,Tposition,beta): 
    '''
    Permet de déterminer si un k-mer peut être étendu. Si c'est possible, renvoie l'indice dans Tposition du k-mer avec lequel il peut être lié, ainsi que ses positions dans les génomes pour lesquels l'extension est possible.
    \param pos : entier, position du k-mer
    \param NKmer : entier, indice du k-mer dans Tposition
    \param Tposition : tableau de tableaux de tableaux d'entiers
    \param beta : entier, seuil de génome à partir duquel on considère la liaison entre k-mers
    \return True si l'extension est possible, False sinon
    \return ksauv : entier, l'indice dans Tposition du k-mer lié
    \return tabRes : tableau de tableaux d'entiers, tableau des positions du k-mer lié dans les différents génomes
    '''
    tabRes = []
    for Kmer in range (0,len(Tposition[0])): #pour chaque k-mer commun du génome de référence
        if pos+1 in Tposition[0][Kmer] : #si une liaison est possible
            ksauv = Kmer
            tabRes.append([pos+1])
            liaison = 1
            genomeVu = []
            #Recherche des autres liaisons possibles
            for genome in range(1,len(Tposition)): #pour chaque génome hors génome de référence
                tabRes.append([])
                for position2 in range (0,len(Tposition[genome][NKmer])):
                    #si une liaison existe dans le génome
                    if (Tposition[genome][NKmer][position2]+1) in Tposition[genome][Kmer]:
                        if genome not in genomeVu :
                            tabRes[-1].append(Tposition[genome][NKmer][position2]+1)
                            liaison += 1
                            genomeVu.append(genome)
            #si le nombre de liaison respecte le seuil beta
            if (liaison>=beta):
                return True,ksauv,tabRes
            else :
                return False,-1,-1
    #s'il n y a pas de liaison possible dans le génome de référence
    return False,-1,-1


def zonesCommunes(Tposition,beta):
    '''
    Crée la structure de données liste_Finale qui contient la succession de sommets constituant les zones communes dans le génome de référence, ainsi que tabCentAdr qui comporte les zones communes pour chaque génome en succession d'arêtes signifiant les liaisons.
    \param Tposition : tableau de tableaux de tableaux d'entiers
    \param beta : entier
    \return liste_Finale : tableau d'entiers
    \return tabCentAdr : tableau de tableaux de tableaux de tableaux d'entiers
    '''
    liste_Finale = []
    listeDesParcourus = []
    tabCentAdr = []
    for i in range (len(Tposition)) :
        tabCentAdr.append([]) 
    res = []
    for i in range (len(Tposition[0])): #pour chaque kmer du génome de référence
        for z in range(len(Tposition[0][i])): #pour chaque position de ce k-mer
            cpt = 1 
            if Tposition[0][i][z] not in listeDesParcourus : #si cette position n'a pas déjè été vue 
                listeDesParcourus.append(Tposition[0][i][z])
                res.append([Tposition[0][i][z]]) #elle devient le début d'une zone commune
                for q in range (len(Tposition)): #on ajoute une zone commune à chaque génome dans tabCentAdr
                    tabCentAdr[q].append([])
                bool,kmer,tabRes = extensible(Tposition[0][i][z],i,Tposition,beta)
                while bool: #tant que l'extension est possible
                    for q in range (len(tabRes)): #pour chaque génome dans tabRes
                        if tabRes[q] != [] :
                            tabCentAdr[q][-1].append([tabRes[q][0]-1,tabRes[q][0]]) #on ajoute à la zone commune qui vient d'être créée pour chaque génome où il y a eu une liaison l'arête de la liaison
                    res[-1].append(tabRes[0][0]) #on ajoute la position étendue à res
                    listeDesParcourus.append(tabRes[0][0])
                    bool,kmer,tabRes = extensible(tabRes[0][0],kmer,Tposition,beta)

    #Réduction de res en ne gardant que les plus grandes zones lorsque qu'il y a chevauchement
    if len(res) > 1 :
        liste_Finale.append(res[0])
        for zone in res[1:]: 
            i = 0
            j = 0
            while i < len(liste_Finale) and zone[j] not in liste_Finale[i]  : #tant qu'on a pas parcouru toutes les zones de liste_Finale ou tant qu'il n y a pas de chevauchement de zone
                j += 1 #on passe à la position suivante
                if j == len(zone) : #si on a parcouru toute la zone 
                    i+=1 #on passe à une autre zone de liste_Finale
                    j = 0 #et on se remet à la première position de zone
            if i == len(liste_Finale) : #si on a pas trouvé de chevauchement, en ayant parcouru toutes les zones de liste_Finale
                liste_Finale.append(zone) #cette zone est unique on l'ajoute à liste_Finale
            else : #si on a trouvé un chevauchement, alors la zone qui a été trouvé est forcément plus complète (*)
                liste_Finale[i] = zone #on écrase la zone précédente avec la nouvelle
    else :
        liste_Finale = res
    liste_Finale.sort()
    return liste_Finale,tabCentAdr
#(*) S'il existe une autre zone où il y a chevauchement non encore traité, alors cette zone est forcément plus grande, c'est-à-dire plus complète. En effet, lorsqu'on crée une zone commune, on vérifie que la première position n'a pas déjà été vue, mais on ne vérifie pas cela lors de l'extension, tout en ajoutant les positions de l'extension comme positions déjà vues.

def redTCA(tabCentAdr):
    '''
    Réduit la structure de données tabCentAdr en éliminant les doublons, ne conservant que les plus grandes zones communes en cas de chevauchement dans le génome de référence.
    \param tabCentAdr : tableau de tableaux de tableaux de tableaux d'entiers, structure qui contient pour chaque génome ses zones communes
    \return tabCentAdrReduit : tableau de tableaux de tableaux de tableaux d'entiers, structure réduite qui contient pour chaque génome ses zones communes
    '''
    tabCentAdrReduit = [] 
    for j in range (len(tabCentAdr)):
        tabCentAdrReduit.append([]) 
    #Même principe que pour la réduction de liste_Finale
    if len(tabCentAdr[0]) > 1 :
        for k in range(len(tabCentAdr)) :
            tabCentAdrReduit[k].append(tabCentAdr[k][0])
        for zIndex,zone in enumerate(tabCentAdr[0]):
            if zone != [] : 
                i = 0
                j = 0
                while i < len(tabCentAdrReduit[0]) and zone[j] not in tabCentAdrReduit[0][i]  :
                    j += 1
                    if j == len(zone) :
                        i+=1 
                        j = 0
                if i == len(tabCentAdrReduit[0]) :
                    for k in range(len(tabCentAdrReduit)) :
                        tabCentAdrReduit[k].append(tabCentAdr[k][zIndex])
                else :
                    for k in range(len(tabCentAdrReduit)) :
                        tabCentAdrReduit[k][i] = tabCentAdr[k][zIndex]
    else :  
        tabCentAdrReduit = tabCentAdr
    return tabCentAdrReduit

def triTCA(tabCentAdr):
    '''
    Trie la structure de données tabCentAdr par ordre croissant selon le génome de référence.
    \param tabCentAdr : tableau de tableaux de tableaux de tableaux d'entiers
    \return tabCentAdrTrie : tableau de tableaux de tableaux de tableaux d'entiers
    '''
    tabCentAdrTrie = []
    for i in range (len(tabCentAdr)) :
        tabCentAdrTrie.append([])
    clas = []
    recup = {}
    #affichMatGenRef(tabCentAdr)
    for groupe_index,groupe in enumerate(tabCentAdr[0]) : #on parcourt les zones communes dans le génome de référence
        #print("\t",groupe)
        recup[groupe[0][0]] = groupe_index #on associe dans le dictionnaire recup la première position de la zone avec son idice de zone
        clas.append(groupe[0][0]) #on ajoute cette position dans le tableau clas
    clas.sort() #on trie par ordre croissant le tableau clas
    for min in clas : 
        for i in range(len(tabCentAdr)) :
            #ajout par ordre croissant dans chaque génome des zones d'indices dont les positions sont les plus petites dans le génome de référence
            tabCentAdrTrie[i].append(tabCentAdr[i][recup[min]]) 
    return tabCentAdrTrie
            
def adr2seq(liste_Finale,saveGkampi):
    '''
    Convertit une liste de positions successives dans le génome de référence en une liste de séquences de nucléotides via le dictionnaire saveGkampi.
    \param liste_Finale : tableau d'entiers
    \param saveGkampi : dictionnaire
    \return listeseq : tableau de chaînes de caractères
    '''
    listeseq = []
    for i in liste_Finale:
        sequence = saveGkampi[i[0]] #Conversion de la première position de la zone en séquence de nucléotides via saveGkampi
        for y in i[1:]:
            sequence = sequence + saveGkampi[y][-1] #Concaténation de la dernière lettre des k-mers aux positions suivantes
        listeseq.append(sequence)
    return listeseq


def sauvGenCons(liste_Finale,listeseq,k,argv):
    '''
    Crée le fichier de sortie fasta qui contient le génome consensus.
    \param liste_Finale : tableau d'entiers, tableau des positions dans le génome de référence des zones communes
    \param listeseq : tableau de chaînes de caractères, tableau des séquences de nucléotides
    \param k : entier, taille du k-mer
    \param argv : tableau de chaînes de caractères, tableau des noms des fichiers
    \return 0
    '''
    scaffoldList = []
    nameliste = []
    infoSup=""
    with open("ResulFasta.fasta","w") as File:
        lecture = SeqIO.to_dict(SeqIO.parse(argv[0],"fasta"))
        #Récupération de tous les scaffolds du génome de référence 
        for seq in lecture:
            if scaffoldList == []:
                scaffoldList.append([0,len(lecture[seq])])
                nameliste.append(seq)
            else :
                scaffoldList.append([(scaffoldList[-1][1])+1,len(lecture[seq])+scaffoldList[-1][1]])
                nameliste.append(seq)
        #Écriture des séquences de listeseq dans le fichier de sortie
        for id,sequence in enumerate(listeseq):
            for compteur,pos in enumerate(scaffoldList):
                infoSup = f" | {liste_Finale[id][0]},{liste_Finale[id][-1] + k-1} VS {pos[0]},{pos[1]} "
                if liste_Finale[id][0]>= pos[0] and liste_Finale[id][-1] + k-1 <= pos[1]: # s'il est contenu dans le scaffold
                    infoSup += "| Pos in seq ("+ nameliste[compteur] +") inside the reference fasta file : begin = " + str( -pos[0] + liste_Finale[id][0])#begin = position de début de la séquence relatif au scaffold
            entete = ">sequence/scaffold_"+ str(id+1) +" | Position with Gkampi index in g1 "+ str(liste_Finale[id][0]) +" ==> "+ str(liste_Finale[id][-1] + k-1) + " | " + str(argv[0]) + infoSup +"\n"
            File.write(entete+sequence+"\n")
    return 0 

def sauvZonesCommunes(liste_Finale):
    '''
    Crée le fichier de sortie resultat.csv qui contient les zones communes pour permettre le débogage.
    \param liste_Finale : tableau d'entiers, tableau des positions dans le génome de référence des zones communes
    \return 0
    '''
    with open("resultat.csv",'w',newline='') as File: 
        Ecriture = csv.writer(File,delimiter=";")
        Ecriture.writerows(liste_Finale)
    return 0 




def affStructSeg(tabCentAdr,tabFrame,app,argv,k,alpha,beta):
    '''
    Crée et affiche la figure matplotlib des zones communes.
    \param tabCentAdr : tableau de tableaux de tableaux de tableaux d'entiers, matrice des zones communes pour chaque génome
    \param tabFrame : objet customtkinter CTkScrollableFrame
    \param app : objet racine de l'interface customtkinter
    \param argv : tableau de chaînes de caractères, noms des fichiers insérés
    \param k : entier, taille du k-mer inséré
    \param alpha : entier, seuil de génome à partir duquel un k-mer commun est considéré
    \param beta : entier, seuil de génome à partir duquel une liaison entre k-mers communs est considérée
    \return 0
    '''
    #tkinter n'accepte que les couleurs en codage hexadécimal or matplotlib donne des couleurs en codage RGB, une fonction de conversion est donc nécessaire
    def rgba_to_hex(rgba):
        '''
        Convertit un codage rgba en un codage hexadécimal.
        \param rgba : codage de couleur en rgba
        \return codage de couleur en hexadécimal
        '''
        r, g, b, _ = rgba
        return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'
    
    def on_closing():
        '''
        Permet de bien fermer la figure matplotlib lorsque l'on ferme l'interface customtkinter.
        '''
        plt.close('all')  #Fermeture de toutes les figures matplotlib
        after_tasks = app.tk.eval("after info").split()
        for after_id in after_tasks: #Fermeture de toutes les fonctionnalités associées à l'objet app
            app.after_cancel(after_id)
        print("fermeture")
        app.destroy()  #Destruction de l'application tkinter

    def tkMatplot(handles,fig):
        '''
        Met en place la fenêtre tkinter qui va héberger la figure matplotlib ainsi que sa légende.
        \param handles : tableau de tuples contenant un codage hexadécimal d'une couleur et une chaîne de caractères, tableau qui contient la légende du schéma
        \param fig : objet de la figure matplotlib
        '''
        def on_click(event,segment_info): 
            '''
            Configure l'affichage de la zone d'un segment lorsqu'il est cliqué.
            \param event : objet de l'événement du clic
            \param segment_info : dictionnaire qui associe à chaque segment sa légende
            '''
            x, y = event.xdata, event.ydata
            if x is not None and y is not None : 
                for line in segment_info : #Itération sur tous les segments du schéma
                    line_xdata = line.get_xdata() 
                    line_ydata = line.get_ydata() 
                    if line_xdata[0] <= x <= line_xdata[1] and line_ydata[0] - 1 <= y <= line_ydata[0] + 1 : #La zone de clic en ordonnée est de 2
                        group_id = segment_info[line] 
                        print(f"Clicked on {group_id}")
                        info_label.config(text=f"Vous avez cliqué sur la zone : {group_id}")
                        break
        root = tabFrame
        def sup() :
            '''
            Configure la suppression de schémas via le bouton associé en fonction du nombre de schémas.
            \param root : objet racine de l'interface customtkinter
            '''
            global nbSub #Récupération de la variable globale nbSub
            root.nametowidget(f"mf{nbSub}").destroy() #Destruction de l'objet Frame qui contient le schéma
            if nbSub == 1 :
                root.nametowidget(f"btn1").destroy() #Si il ne reste que le premier schéma : destruction du bouton de suppression
            nbSub -= 1 
        
        global nbSub 
        print(f"Nombre de schéma : {nbSub}")
        main_frame = ttk.Frame(root, padding='0.05i', name=f"mf{nbSub}")#Création de l'objet Frame qui va contenir le schéma

        #Création des labels d'information qui renseignent les paramètres du schéma affiché
        name = "Reference genome = "+argv[0]+ "| compared genomes = "
        for fic in argv[1:] :
            name += fic + " "        
        labelF = ttk.Label(main_frame, text= f"{name} ", background='white')
        labelF2 = ttk.Label(main_frame,text=f"| k-mer size = {k} | genome with the same k-mer = {alpha} | genome with the same link between k-mer = {beta}", background='white')
        #Placement des labels créés dans main_frame
        labelF.pack(fill=tk.BOTH, expand=1)
        labelF2.pack(fill=tk.BOTH, expand=1)
 
        
        #Création de l'étiquette qui va annoncer le clic sur segment
        info_label = tk.Label(main_frame, text="Cliquez sur un segment pour voir les informations du groupe", bg="white", fg="black", font=("Arial",12))
        info_label.pack(side=tk.BOTTOM, fill=tk.X) #Cette étiquette est placée en bas de tout et prend toute la largeur

        #Création d'un objet Canvas qui va contenir la figure matplotlib
        scroll_canvas = tk.Canvas(main_frame, background="blue")
        #Création des objets Scrollbar qui vont permettre son défilement
        scrollbarH = ttk.Scrollbar(main_frame, orient=tk.HORIZONTAL, command=scroll_canvas.xview)
        scrollbarH.pack(side=tk.BOTTOM, fill=tk.X)
        scrollbarV = ttk.Scrollbar(main_frame, orient=tk.VERTICAL, command=scroll_canvas.yview)
        scrollbarV.pack(side=tk.RIGHT, fill=tk.Y)
        #Placement de main_frame et scroll_canvas
        main_frame.pack(side=tk.TOP, expand=1, fill=tk.BOTH)
        scroll_canvas.pack(side=tk.LEFT,  expand=1, fill=tk.BOTH)
        #Configuration du défilement de scroll_canvas
        scroll_canvas.configure(xscrollcommand=scrollbarH.set,yscrollcommand=scrollbarV.set)
        scrollbarH.config(command=scroll_canvas.yview)
        scrollbarV.config(command=scroll_canvas.xview)
        #Création de canvas_frame qui va contenir la figure matplotlib
        canvas_frame = ttk.Frame(scroll_canvas)
        scroll_canvas.create_window((0,0), window=canvas_frame,anchor='center', height=300)
        canvas_frame.pack(expand=1,fill=tk.BOTH)
        canvas = FigureCanvasTkAgg(fig, master=canvas_frame)
        canvas.draw()
        canvas_widget = canvas.get_tk_widget()
        canvas_widget.pack(side=tk.LEFT, fill=tk.BOTH, expand=1) 
        canvas_frame.update_idletasks()#Mise à jours de canvas_frame 
        scroll_canvas.config(scrollregion=scroll_canvas.bbox("all"))#Configuration du défilement sur tout le contenu de scroll_canvas

        #Ajout de la Toolbar de matplotlib en widget tkinter
        toolbar = NavigationToolbar2Tk(canvas, main_frame) 

        #Création de l'objet Frame de la légende du schéma
        legend_frame = ttk.Frame(main_frame) 
        legend_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=1) 
        scrollbar = ttk.Scrollbar(legend_frame, orient=tk.VERTICAL) 
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        legend_canvas = tk.Canvas(legend_frame) 
        legend_frame_inner = ttk.Frame(legend_canvas)
        scrollbar.config(command=legend_canvas.yview)
        legend_canvas.create_window((0, 0), window=legend_frame_inner, anchor='nw') 
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

        
        #Mise à jours de legend_frame et scroll_canvas avec configuration du défilement sur tout leur contenu
        legend_frame_inner.update_idletasks() 
        legend_canvas.config(scrollregion=legend_canvas.bbox("all"))     
        scroll_canvas.update_idletasks()
        scroll_canvas.config(scrollregion=scroll_canvas.bbox("all")) 

        canvas.mpl_connect("button_press_event",on_click) 
        if nbSub == 1 :
            #Création du bouton de suppression de schémas au premier submit
            btn1 = ttk.Button(root, text ="Supprimer le dernier schéma", command = sup, name= f"btn{nbSub}")
            btn1.pack(pady = 10) 

    #-------Création de la figure matplotlib-----------
    #Configuration des couleurs associés aux zones communes via la bibliothèque "viridis"
    num_colors = max(len(data_set) for data_set in tabCentAdr) 
    cmap = plt.cm.get_cmap('viridis', num_colors)
    #Récupération de la plus petite et de la plus grande valeur pour l'échelle du schéma
    min_val = min(min(zone[0] for group in data_set for zone in group) for data_set in tabCentAdr)
    max_val = max(max(zone[1] for group in data_set for zone in group) for data_set in tabCentAdr)
    #Mise en place de l'écart minimal pour chaque génome
    max_data_set = max( len(data_set) for data_set in tabCentAdr)
    #Création du schéma
    fig, ax = plt.subplots(figsize=(5,2))
    legend_handles = []
    segment_info = {}
    y_mid = []
    labels = []
    for ds_index, data_set in enumerate(tabCentAdr):
        y_value =3*max_data_set*ds_index #L'ordonnée de début pour placer les zones communes de chaque génome
        y_mid.append(y_value+(len(data_set)/2)) #L'ordonnée à laquelle on place la modalité (Génome i)
        for group_index, group in enumerate(data_set):
            y_value += 1 #Écart pour éviter tout chevauchement de zone
            color = cmap(group_index/num_colors) 
            hex_color = rgba_to_hex(color) 
            if ds_index == 0 : 
                #Récupère les bornes inférieure et supérieure du groupe dans le génome de référence, ce qui va servir pour la légende.
                min_group = group[0][0] 
                max_group = group[-1][1]
                labels.append(f'{min_group}-{max_group}')
                legend_handles.append((hex_color,f'{min_group}-{max_group+k-1}')) 
            for zone in group:
                #Traçage des zones du groupe
                line, = ax.plot(zone, [y_value,y_value], linewidth=10, color=color) 
                segment_info[line] = labels[group_index] #Association de chaque segment avec sa légende via le dictionnaire segment_info
    app.protocol("WM_DELETE_WINDOW", on_closing)
    
    #Mise en place des modalités
    ax.set_xlim(min_val - 1, max_val + 1) 
    ax.set_yticks(y_mid)
    ystickslables = []
    for fic_index,fic in enumerate(argv):
        ystickslables.append(f"Génome {fic_index+1} in {fic}") 
    ax.set_yticklabels(ystickslables)
    fig.tight_layout(rect=[0, 0, 1, 1])
    tkMatplot(legend_handles,fig)
    return 0
    


def main(argv,kmer,alpha,beta,tabFrame2,app) :
    '''
    La fonction qui lance le programme.
    \param argv : tableau de chaînes de caractères 
    \param kmer : chaîne de caractères
    \param alpha : chaîne de caractères
    \param beta : chaîne de caractères
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
    affichMatGenCompMeilleur(triTCA(redTCA(tabCentAdr)))
    listeSeq = adr2seq(liste_Finale,saveGKampi)

    #Les sorties finales
    sauvGenCons(liste_Finale,listeSeq,kmer,argv) 
    sauvZonesCommunes(liste_Finale)

    #Affichage du schéma des zones communes
    affStructSeg(triTCA(redTCA(tabCentAdr)),tabFrame2,app,argv,kmer,alpha,beta)

