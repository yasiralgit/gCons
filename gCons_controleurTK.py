from gCons_fonctionsTK import *
from customtkinter import * 

root = CTk()
root.geometry("500x400")

set_default_color_theme("green")

tabview = CTkTabview(master=root)
tabview.pack(expand= 1, anchor='n',padx=20,pady=20,fill="both")

tabview.add("Tab 1")
tabview.add("Tab 2")


label = CTkLabel(master=root, text="Zones Communes")
label.place(relx=0.5, anchor="n")

frame = CTkScrollableFrame(master=tabview.tab("Tab 1"), border_width=2, orientation="vertical")
frame.pack(expand=True, anchor="n",pady=30,fill='both')

frameT2 = CTkScrollableFrame(master=tabview.tab("Tab 2"), border_width=2, orientation="vertical")
frameT2.pack(expand=True, anchor="n",pady=30, fill='both')

frameInnit = CTkFrame(master=frame, border_width=2)
frameInnit.pack(expand=True, anchor="n",pady=30, padx=10)

frame2 = CTkFrame(master=frame, border_width=2)
frame2.pack(expand=True, anchor="n",pady=30, padx=10)

frame3 = CTkFrame(master=frame, border_width=2)
frame3.pack(expand=True, anchor="n",pady=30,padx=10)

frame4 = CTkFrame(master=frame, border_width=2)
frame4.pack(expand=True, anchor="n",pady=30,padx=10)




def click_hander(): 
    '''
    Lance le programme avec les informations renseignées une fois le bouton submit appuyé
    '''
    print(f"Fichier(s) inséré(s) : {textbox.get()}")
    print(f"Taille de k-mer souhaitée : {entry1.get()}")
    print(f"Pourcentage alpha : {entryA.get()}")
    print(f"Pourcentage beta : {entryB.get()}")
    argv = textbox.get().split()
    main(argv, entry1.get(),entryA.get(),entryB.get(),frameT2,root)
global textbox
global entry1
global entryA
global entryB 

label = CTkLabel(master=frameInnit, text="Saisissez le nom des fichiers fasta contenant les génomes : ")
textbox = CTkEntry(master=frameInnit,width=300, text_color="#FFCC70")

label.pack(anchor="w", expand=True, padx= 10, pady = 10)
textbox.pack(anchor="center", expand=True, pady=10, padx= 30)



label = CTkLabel(master=frame2, text="Saisissez la taille de k-mer souhaitée")
entry1 = CTkEntry(master=frame2, width=50, text_color="#FFCC70")
btn = CTkButton(master=frame,text="Submit",command=click_hander)

label.pack(anchor="w", expand=True, padx= 10, pady = 10)
entry1.pack(anchor="center", expand=True, pady=10, padx= 30)
btn.pack(anchor="s",expand=True, padx=30, pady= 20)



alpha = CTkLabel(master=frame3, text="Saisissez le pourcentage de génome ayant le meme k-mer")
entryA = CTkEntry(master=frame3, width=50, text_color="#FFCC70")

alpha.pack(anchor="w", expand=True, padx= 10, pady = 10)
entryA.pack(anchor="center", expand=True, pady=10, padx= 30)

beta = CTkLabel(master=frame4, text="Saisissez le pourcentage de génome ayant la meme liaison de k-mer")
entryB = CTkEntry(master=frame4, width=50, text_color="#FFCC70")

beta.pack(anchor="w", expand=True, padx= 10, pady = 10)
entryB.pack(anchor="center", expand=True, pady=10, padx= 30)




root.mainloop()





