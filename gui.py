from tkinter import *
from tkinter import filedialog

from toddCox import *
from tkinter.scrolledtext import ScrolledText
from tkinter import ttk #treeview

class MyWindow:
    def __init__(self, win):
        self.table = []
        self.x = -1
        self.lbGroup=Label(win, text='Group G=<S,R>')
        self.lbGroup.place(x=50, y=20)
        
        #####
        self.lbl1=Label(win, text='Generators:')
        self.lbl1.place(x=60, y=50)
        
        self.lbl2=Label(win, text='Relators:')
        self.lbl2.place(x=65, y=75)
        
        self.lbH=Label(win, text='Group H=<Y>')
        self.lbH.place(x=50, y=120)
        
        self.lbl3=Label(win, text='Generators')
        self.lbl3.place(x=60, y=150)
        ####
        

        #Entrada de arriba, derecha de generators G
        self.t1=Entry(win, bd=2)
        self.t1.place(x=140, y=50)
        
        #Entrada media, derecha de relators G
        self.t2=Entry(win,bd=2)        
        self.t2.place(x=140, y=75)

        #Entrada de abajo, derecha de generators H
        self.t3=Entry(win, bd=2)
        self.t3.place(x=140, y=150)
        
        #index
        self.ti=Entry(win, bd=2)
        self.ti.place(x=420, y=20)

        #G
        self.tg=Entry(win, bd=2)
        self.tg.place(x=420, y=50)
        
        #H
        self.th=Entry(win, bd=2)
        self.th.place(x=420, y=80)
    
        self.lbGroup=Label(win, text='You can browse a file instead:')
        self.lbGroup.place(x=50, y=200)
        
        self.explore = Button(win, text = "Browse File", bd=2, command = self.browseFiles) 
        self.explore.place(x=85, y=226)
      
        self.label_file_explorer = Label(win, text = "File Explorer", fg = "blue") 
        self.label_file_explorer.place(x=50, y=260)

        
        self.validate = Button(win, text = "Apply Todd Coxeter", bd=2,
                              fg="green", command = self.validate) 
        self.validate.place(x=200, y=226)
        
        
        self.lbGroup=Label(win, text='Index of H in G: ')
        self.lbGroup.place(x=320, y=20)
        
        self.lbGroup=Label(win, text='Order of G: ')
        self.lbGroup.place(x=329, y=50)
        
        self.lbGroup=Label(win, text='Order of H: ')
        self.lbGroup.place(x=329, y=80)
        
        
        
        self.lbGroup=Label(win, text='Cosets: ')
        self.lbGroup.place(x=340, y=130)
     
        
        
        self.tick = Checkbutton(win, text="Show Cosets", variable=IntVar(), command=self.hello)
        self.tick.place(x=495, y=130)
        
        
        
        #tabla de cosets
        

        self.xscrollbar = Scrollbar(win, orient=HORIZONTAL)
        self.xscrollbar.place(x=534, y=323)

        self.yscrollbar = Scrollbar(win)
        #self.yscrollbar.grid(row=3, column=4, sticky=N+S+E+W)

        self.text = ScrolledText(win, wrap=NONE, width=30, height=10,
                    xscrollcommand=self.xscrollbar.set,
                    yscrollcommand=self.yscrollbar.set)
        #self.text.grid(row=17, column=16)
        self.text.place(x=340, y=160)
        #self.text.configure(state="disabled")

        self.xscrollbar.config(command=self.text.xview)
        self.yscrollbar.config(command=self.text.yview)
        
        
        
        #self.scrollH = Scrollbar(win, orient='horizontal')
        #self.text = ScrolledText(win, width=30, height=10, borderwidth=3,
                                 #relief="sunken", xscrollcommand=self.scrollH.set)
        #self.text.place(x=340, y=160)
        #self.scrollH.pack(side="bottom", fill="x")
        #self.scrollH.config(command=self.text.xview)
        
        
    def hello(self):
        self.x=self.x+1
        if self.x%2==0:
            print(len(self.table))
            if len(self.table)<=25:
                self.text.delete("1.0",'end')
                self.text.insert(END, self.table)
            else:
                self.text.insert(END, "table of cosets is too big")
        else:
            self.text.delete("1.0",'end')
        
    def validate(self):
        t1= self.t1.get()
        t2= self.t2.get()
        t3 = self.t3.get()
        
        rep = {",": " ", "=": " ", "1": " ", "{": " ", "}": " "}

        gen = replace_all(t1, rep).split()
        rels = replace_all(t2, rep).split()
        genH = replace_all(t3, rep).split()
        
        
        #calcular índice
        nueva = CosetTable(2*len(gen), gen, rels, genH)
        nueva.CosetEnumeration()
        n = nueva.finalCosets()
        print(n)
        
        self.ti.delete(0, 'end')
        self.ti.insert(END, str(n))
        
        
        
        #calcular |G|
        if len(genH) == 0:
            g = n
        else:
            group = CosetTable(len(gen)*2, gen, rels, [])
            group.CosetEnumeration()
            g = group.finalCosets()
        
        self.tg.delete(0, 'end')
        self.tg.insert(END, str(g))
        
        
        h = int(g/n)
        self.th.delete(0, 'end')
        self.th.insert(END, str(h))
        
        
        #table of classes
        #nueva.pretty_print()
        self.table = list(range(0,n))
        
        if n>25:
            print("cant")
            self.text.delete("1.0",'end')
            if self.x%2==0:
                self.text.insert(END, "table of cosets is too big")
        else:
            nueva.pretty_print()
            self.table = nueva.table
            self.text.delete("1.0",'end')
            
            if self.x%2==0:
                self.text.insert(END, self.table)


        
        
        
    def browseFiles(self): 
        filename = filedialog.askopenfilename(initialdir = "/", 
                    title = "Select a File",filetypes = (("Text files", 
                                        "*.txt*"), ("all files","*.*"))) 
        content = open(filename, 'r')
        # Change label contents 
        self.label_file_explorer.configure(text="File Opened: "+filename) 

       
        gen = content.readline()
        rels = content.readline()
        genH = content.readline()
        
        self.t1.delete(0, 'end')
        self.t2.delete(0, 'end')
        self.t3.delete(0, 'end')
        
        self.tg.delete(0, 'end')
        self.th.delete(0, 'end')
        self.ti.delete(0, 'end')
        self.text.delete("1.0",'end')
        
        self.t1.insert(END, gen)
        self.t2.insert(END, rels)
        self.t3.insert(END, genH)
        
       

window=Tk()
mywin=MyWindow(window)


window.title('Todd Coxeter Algorithm')
window.geometry("650x380+300+100")
#window.configure(background="white")
window.mainloop()