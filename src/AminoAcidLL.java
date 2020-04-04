class AminoAcidLL {
  char aminoAcid;
  String[] codons;
  int[] counts;
  AminoAcidLL next;

  AminoAcidLL() {

  }


  /********************************************************************************************/
  /* Creates a new node, with a given amino acid/codon
   * pair and increments the codon counter for that codon.
   * NOTE: Does not check for repeats!! */
  AminoAcidLL(String inCodon) {
    this.aminoAcid = AminoAcidResources.getAminoAcidFromCodon(inCodon);
    this.codons = AminoAcidResources.getCodonListForAminoAcid(aminoAcid);
    this.next = null;
  }

  private void incrCodons(String inCodon) {
    for (int i = 0; i < codons.length; i++) {
      if (codons[i].equals((inCodon))){
      counts[i]++;
    }
  }
}
  /********************************************************************************************/
  /* Recursive method that increments the count for a specific codon:
   * If it should be at this node, increments it and stops, 
   * if not passes the task to the next node. 
   * If there is no next node, add a new node to the list that would contain the codon. 
   */
    private void addCodon (String inCodon) {
        if (aminoAcid == AminoAcidResources.getAminoAcidFromCodon(inCodon)) {
          incrCodons(inCodon);
        }
        else if (next != null) {
          next.addCodon(inCodon);
        }
        else{
          next = new AminoAcidLL(inCodon);
          addCodon(inCodon);
      }
    }

  /********************************************************************************************/
  /* Shortcut to find the total number of instances of this amino acid */
  private int totalCount(){
    int sum = 0;
    for(int i = 0; i < counts.length; i++){
      sum += counts[i];
    }
    return  sum;
  }

  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int totalDiff(AminoAcidLL inList){
    return Math.abs(totalCount() - inList.totalCount());
  }


  /********************************************************************************************/
  /* helper method for finding the list difference on two matching nodes
  *  must be matching, but this is not tracked */
  private int codonDiff(AminoAcidLL inList){
    int diff = 0;
    for(int i=0; i<codons.length; i++){
      diff += Math.abs(counts[i] - inList.counts[i]);
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursive method that finds the differences in **Amino Acid** counts. 
   * the list *must* be sorted to use this method */
  public int aminoAcidCompare(AminoAcidLL inList){
    if(!inList.isSorted()){
      return 0;
    }
    int diff = 0;
    if(inList == null){
      diff += totalCount();
    }
    if(next != null){
      diff += next.aminoAcidCompare(inList);
    }
    else if(aminoAcid == inList.aminoAcid){
      diff += totalDiff(inList);
    }
    if(next != null){
      diff += aminoAcidCompare(inList.next);
    }
    if(next == null && inList.next != null){
      diff += aminoAcidCompare(inList.next);
    }
    else if(next != null && aminoAcid < inList.aminoAcid){
      diff += totalCount();
      if(next != null){
        diff += next.aminoAcidCompare(inList);
      }
      else if(next == null || aminoAcid > inList.aminoAcid){
        diff += inList.totalCount();
        if(inList.next != null){
          diff += aminoAcidCompare(inList.next);
        }
      }
    }
    return diff;
  }

  /********************************************************************************************/
  /* Same ad above, but counts the codon usage differences
   * Must be sorted. */
  public int codonCompare(AminoAcidLL inList) {
    if (inList.isSorted()) {
      return 0;
    }
    int diff = 0;
    if (inList == null) {
      diff += totalCount();
      next.codonCompare(inList);
    }
    if (next != null) {
      diff += next.codonCompare(inList.next);
    } else if (aminoAcid == inList.aminoAcid) {
        diff += codonDiff(inList);

      if (next != null) {
        diff += next.codonCompare(inList.next);
      }
      if(next == null && inList.next != null){
        diff += codonCompare(inList.next);
      }
      else if(next != null && aminoAcid < inList.aminoAcid){
        diff += totalCount();
        if(next != null){
          diff += next.codonCompare(inList);
        }
        else if(next == null || aminoAcid > inList.aminoAcid){
          diff += inList.totalCount();
          if(inList.next != null){
            diff += codonCompare(inList.next);
          }
        }
      }
    }
    return diff;
  }

  /********************************************************************************************/
  /* Recursively returns the total list of amino acids in the order that they are in in the linked list. */
  public char[] aminoAcidList() {
    if (next == null) {
      return new char[]{aminoAcid};
    }
    char[] a = next.aminoAcidList();
    char[] finalA = new char[a.length + 1];
    finalA[0] = aminoAcid;

    for (int i = 1; i < finalA.length; i++) {
      finalA[i] = a[i - 1];
    }
    return finalA;
  }

  /********************************************************************************************/
  /* Recursively returns the total counts of amino acids in the order that they are in in the linked list. */
  public int[] aminoAcidCounts(){
    if(next == null){
      return new int[]{};
    }
    else{
      int[] node = next.aminoAcidCounts();
      int[] counts = new int[node.length + 1];

      for(int i = 0; i < node.length; i++){
        counts[0] = aminoAcid;
        counts[i + 1] = node[i];
      }
      return aminoAcidCounts();
    }
  }


  /********************************************************************************************/
  /* recursively determines if a linked list is sorted or not */
  public boolean isSorted(){
    if(next == null){
      return true;
    }
    if(next.aminoAcid > aminoAcid){
      return false;
    }
    return next.isSorted();
  }


  /********************************************************************************************/
  /* Static method for generating a linked list from an RNA sequence */
  public static AminoAcidLL createFromRNASequence(String inSequence){
    AminoAcidLL block = new AminoAcidLL(inSequence.substring(0, 3));
    boolean test = true;

    if(inSequence.substring(0, 3).charAt(0) == '*'){
      block.addCodon(inSequence.substring(0, 3));
      test = false;
    }
    else{
      block.addCodon(inSequence.substring(0, 3));
    }

    for (int i = 3; i < inSequence.length() - 2; i += 3) {
      //use string.subString(start, end)
      if (inSequence.charAt(i) == '*') {
        block.addCodon(inSequence.substring(i, i + 3));
        test = false;
      }
      else{
        block.addCodon(inSequence.substring(i, i + 3));
      }
    }
    return block;
  }


  /********************************************************************************************/
  /* sorts a list by amino acid character*/
  public static AminoAcidLL sort(AminoAcidLL inList){
    if(inList.isSorted()){
      return inList;
    }
    else{
      for(AminoAcidLL i = inList; i.next != null; i = i.next){
        for(AminoAcidLL j = i.next; j.next != null; j = j.next){
          if (i.aminoAcid > j.aminoAcid) {
            AminoAcidLL temperature = i;
            j.next = i;
            i = temperature;
          }
        }
      }
    }
    return null;
  }
}