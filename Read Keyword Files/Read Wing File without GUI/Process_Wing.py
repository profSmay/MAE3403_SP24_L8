from Wing_Class import Wing

def GetWing():
    # get the filename using the OPEN dialog
    filename = 'Wing Input File - good Loads.txt'
    #filename = 'Wing Input File - bad Loads.txt'
    f1 = open(filename, 'r')  # open the file for reading
    data = f1.readlines()  # read the entire file as a list of strings
    f1.close()  # close the file  ... very important

    wing = Wing()  # create a wing instance (object)

    try:
        wing.processWingData(data)
        return wing
    except:
        print('Bad File')

def main():
    wing=GetWing()
    print(wing.title)
    print(wing.sparH)
    print(wing.sparW)
    print(wing.up_dist)
    print(wing.aft_dist)
    print(wing.up_point)
    print(wing.aft_point)

main()