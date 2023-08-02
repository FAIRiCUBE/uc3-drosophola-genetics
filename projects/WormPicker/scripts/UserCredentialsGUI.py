from tkinter import Tk, Label,Entry, IntVar, Checkbutton, Button
import tkinter.messagebox as messagebox

username_1=""
password_1=""
hostadress_1=""
default_host="https://fairicube.rasdaman.com/rasdaman/ows"

def get_input():
    username = entry_username.get()
    global username_1
    username_1 = username
    global password_1
    password = entry_password.get()
    password_1 = password
    global hostadress_1
    hostadress = entry_host.get()
    hostadress_1=hostadress
    window.destroy()

window = Tk()
window.title('User Credentials')

label_username = Label(window, text='Username:')
label_username.pack(padx=15, pady=5)

entry_username = Entry(window)
entry_username.pack(padx=15, pady=5)

label_password = Label(window, text='Password:')
label_password.pack(padx=15, pady=6)

entry_password = Entry(window, show='*')  # Use show='*' to hide password input
entry_password.pack(padx=15, pady=7)

label_host= Label(window, text="Host:")
label_host.pack(padx=15, pady=5)

entry_host= Entry(window)
entry_host.insert(0, default_host)
entry_host.pack(padx=15, pady=7)

btn_submit = Button(window, text='Submit', command=get_input)
btn_submit.pack(pady=10)
window.mainloop()

#RASDAMAN_CRED_USERNAME=username_1
#RASDAMAN_CRED_PASSWORD=password_1 

print("Username:", username_1)
print("Password:", password_1)

with open("/media/ssteindl/fairicube/uc3/uc3-drosophola-genetics/.env3", "w") as file:
    # Convert variables to strings and write them to the file
    file.write("RASDAMAN_CRED_USERNAME=" + username_1 + "\n")
    file.write("RASDAMAN_CRED_PASSWORD=" + password_1 + "\n")
    file.write("RASDAMAN_SERVICE_ENDPOINT=" + hostadress_1 + "\n")


def show_selected_options():
    selected_option1 = 0
    selected_option2 = 0
    selected_option3 = 0
    if var_option_a.get() == 1:
        selected_option1 = layer1
    if var_option_b.get() == 1:
        selected_option2 = layer2
    if var_option_c.get() == 1:
        selected_option3 = layer3
    options=[selected_option1, selected_option2, selected_option3]
    selection=[element for element in options if element != 0]
    if any(options):
        for x in (selection):
            messagebox.showinfo("Selected Option", f"You chose {x}")
    else:
        messagebox.showinfo("Selected Option", "Please select an option.")

window = Tk()
window.title("Select an Option")

layer1="TEMPERATURE"
layer2="HUMIDITY"
layer3="IMPERVIOUSNESS"

label_option = Label(window, text="Select an option:")
label_option.pack()

var_option_a = IntVar()
check_option_a = Checkbutton(window, text=layer1, variable=var_option_a)
check_option_a.pack()

var_option_b = IntVar()
check_option_b = Checkbutton(window, text=layer2, variable=var_option_b)
check_option_b.pack()

var_option_c = IntVar()
check_option_c = Checkbutton(window, text=layer3, variable=var_option_c)
check_option_c.pack()

btn_finished = Button(window, text="Finished", command=show_selected_options)
btn_finished.pack(pady=10)

window.mainloop()


