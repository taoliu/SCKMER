#!/usr/bin/env python3

def is_simplerepeat ( kmer ):
    l = len(kmer)
    for simple_repeat_l in range( 1, l//2 + 1 ):
        is_repetitive = True
        repeat = kmer[:simple_repeat_l]
        for check_s in range( simple_repeat_l, l - l % simple_repeat_l, simple_repeat_l ):
            chunk = kmer[ check_s: check_s+simple_repeat_l ]
            is_repetitive = is_repetitive and (repeat == kmer[ check_s: check_s+simple_repeat_l ])
        # check the last chunk
        chunk = kmer[ check_s + simple_repeat_l: ]
        is_repetitive = is_repetitive and ( repeat[ : len(chunk) ] == chunk )
        if is_repetitive:
            #print (f"{kmer} simple repeat {repeat}")
            return True
    return False


def is_too_simple ( kmer, max_diff_bp = 2 ):
    nA = kmer.count("A")
    nC = kmer.count("C")
    nT = kmer.count("T")
    nG = kmer.count("G")
    l = len( kmer ) - max_diff_bp
    if nA >= l or nC >= l or nT >= l or nG >= l:
        return True
    else:
        return False
    
def has_same_bp ( kmer, min_length = 8 ):
    A = "A" * min_length
    C = "C" * min_length
    T = "T" * min_length
    G = "G" * min_length
    if kmer.find(A) != -1 or kmer.find(C) != -1 or kmer.find(T) != -1 or kmer.find(G) != -1:
        return True
    else:
        return False
    

def main():
    with open("all.uniq.kmer.txt", "r") as fhd:
    #with open("test.txt", "r") as fhd:
        l = fhd.readline()
        print (l.rstip())
        for l in fhd:
            kmer = l.rstrip()
            if is_too_simple ( kmer, max_diff_bp = 2 ):
                #print (f"{kmer} is too simple!")
                continue 
            elif is_simplerepeat ( kmer ):
                #print (f"{kmer} contains simple repeats!")
                continue
            elif has_same_bp ( kmer ):
                #print (f"{kmer} contains same bps!")
                continue            
            else:
                #print (f"{kmer} pass test!")
                print( kmer )
    



        
