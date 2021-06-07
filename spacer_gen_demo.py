from hetero_spacer_generator import HeteroGen
from Bio.Seq import Seq

valid_input = False
inpt = ''
while not valid_input:
    inpt = input("Please enter you sequence: ")

    if len(inpt) < 0 or inpt is None:
        print("Invalid Input...\n")
        continue

    for char in inpt:
        if char not in "ATCG":
            print("Invalid Input...\n")
            continue
    valid_input = True

seq = Seq(inpt)

valid_input = False
while not valid_input:
    inpt = input("Please enter the max heterogeneity spacer length: ")
    if not inpt.isnumeric():
        print("Invalid Input...\n")
        continue
    if int(inpt) > len(seq):
        print("Spacer cannot be longer than sequence...\n")
        continue
    valid_input = True

print("Will now generate spacer regions ensuring 12 bases of heterogeneity. ")

hg = HeteroGen(max_spacer_length=int(inpt))

spacer_combos = hg.get_all_spacer_combos(seq)

if spacer_combos == []:
    print("Sorry! No spacer combos were found that match your specifications.")
    exit(0)

inpt = input(
    ''.join(["Found ", str(len(spacer_combos)), " combinations of spacers "
                                                "matching "
                                           " your specifications. Enter "
                                           "\"View\" to view."]))
if inpt == "View" or inpt == "view":
    print("\"+\" denotes a base within the spacer region.")
    hg.visualise_spacer_combos(spacer_combos, seq)

while True:
    valid_input = False
    inpt = ''
    while not valid_input:
        inpt = input("Please enter you sequence: ")

        if len(inpt) < 0 or inpt is None:
            print("Invalid Input...\n")
            continue

        for char in inpt:
            if char not in "ATCG":
                print("Invalid Input...\n")
                continue
        valid_input = True

    seq = Seq(inpt)

    valid_input = False
    while not valid_input:
        inpt = input("Please enter the max heterogeneity spacer length: ")
        if not inpt.isnumeric():
            print("Invalid Input...\n")
            continue
        if int(inpt) > len(seq):
            print("Spacer cannot be longer than sequence...\n")
            continue
        valid_input = True

    print("Will now generate spacer regions ensuring 12 bases of heterogeneity. ")

    hg = HeteroGen(max_spacer_length=int(inpt))

    spacer_combos = hg.get_all_spacer_combos(seq)

    if spacer_combos == []:
        print("Sorry! No spacer combos were found that match your specifications.")
        exit(0)

    inpt = input(
        ''.join(["Found ", str(len(spacer_combos)), " combinations of spacers "
                                                    "matching "
                                                    " your specifications. Enter "
                                                    "\"View\" to view."]))
    if inpt == "View" or inpt == "view":
        print("\"+\" denotes a base within the spacer region.")
        hg.visualise_spacer_combos(spacer_combos, seq)

