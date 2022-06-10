MODULE Poseidon_Letters_Module



IMPLICIT NONE

CHARACTER(LEN=1), DIMENSION(1:26), PUBLIC, PARAMETER  ::                &
    Letter_Table_Upper = (/ "A","B","C","D","E","F","G","H","I","J",    &
                            "K","L","M","N","O","P","Q","R","S","T",    &
                            "U","V","W","X","Y","Z" /)

CHARACTER(LEN=1), DIMENSION(1:26), PUBLIC, PARAMETER  ::                &
    Letter_Table_Lower = (/ "a","b","c","d","e","f","g","h","i","j",    &
                            "k","l","m","n","o","p","q","r","s","t",    &
                            "u","v","w","x","y","z" /)

END MODULE  Poseidon_Letters_Module

