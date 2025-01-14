from tkinter import Tk, Label,Entry, IntVar, Checkbutton, Button
import tkinter.messagebox as messagebox

#username_1=""
#password_1=""
#hostadress_1=""
#default_host="https://fairicube.rasdaman.com/rasdaman/ows"

def get_input(entry_username, entry_password, entry_host, window,path):
    username = entry_username.get()
    global username_1
    username_1 = username
    global password_1
    password = entry_password.get()
    password_1 = password
    global hostadress_1
    hostadress = entry_host.get()
    hostadress_1=hostadress
    with open(path, "w") as file:
    # Convert variables to strings and write them to the file
        file.write("RASDAMAN_CRED_USERNAME=" + username_1 + "\n")
        file.write("RASDAMAN_CRED_PASSWORD=" + password_1 + "\n")
        file.write("RASDAMAN_SERVICE_ENDPOINT=" + hostadress_1 + "\n")
    window.destroy()

def saveCredentials(path):
    window = Tk()
    window.title('User Credentials for HOST')
    default_host="https://fairicube.rasdaman.com/rasdaman/ows"
    label_username = Label(window, text='Username:')
    label_username.pack(padx=15, pady=5)
    entry_username = Entry(window)
    entry_username.pack(padx=15, pady=5)
    label_password = Label(window, text='Password:')
    label_password.pack(padx=15, pady=6)
    entry_password = Entry(window, show='*')  # Use show='*' to hide password input
    entry_password.pack(padx=15, pady=7)
    label_host= Label(window, text="HostAdress:")
    label_host.pack(padx=15, pady=5)
    entry_host= Entry(window)
    entry_host.insert(0, default_host)
    entry_host.pack(padx=15, pady=7)
    btn_submit = Button(window, text='Submit', command=lambda: get_input(entry_username, entry_password, entry_host, window,path))
    btn_submit.pack(pady=10)
    window.mainloop()


