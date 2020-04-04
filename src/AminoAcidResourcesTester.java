import static org.junit.jupiter.api.Assertions.*;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;
class AminoAcidResourcesTester{
  @Test
  public void allCodons(){
    System.out.print("hello");
    char[] rna = {'A','C','U','G'};
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    for(int i=0; i<4; i++){
      for(int j=0; j<4; j++){
        for(int k=0; k<4;k++){
          String s = new String(new char[]{rna[i],rna[j],rna[k]});
          char aaOut = AminoAcidResources.getAminoAcidFromCodon(s);
          if(aaOut != '*'){
            String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aaOut);
            boolean found = false;
            for(int l=0; l<codonList.length; l++){
              found |= (codonList[l].equals(s));
            }
            if(!found) System.err.println("Codon " + s + " not found, said AA was " + aaOut);
          }

          aaOut = AminoAcidResources.getAminoAcidFromCodon(s.toLowerCase());
          if(aaOut != '*'){
            String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aaOut);
            boolean found = false;
            for(int l=0; l<codonList.length; l++){
              found |= (codonList[l].equals(s));
            }
            if(!found) System.err.println("Codon " + s + " not found, said AA was " + aaOut);
          }
        }
      }
    }
  }
  @Test
  public void allAAs(){
    System.out.print("hello");
    char[] aa = {'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W'};
    for(int i=0; i<aa.length; i++){
      String[] codonList = AminoAcidResources.getCodonListForAminoAcid(aa[i]);
      for(int l=0; l<codonList.length; l++){
        if(aa[i] != AminoAcidResources.getAminoAcidFromCodon(codonList[l])){
          System.err.println("AA " + aa[i] + " not found, said codon was " + codonList[l]);
        }
      }
      codonList = AminoAcidResources.getCodonListForAminoAcid(Character.toLowerCase(aa[i]));
      for(int l=0; l<codonList.length; l++){
        if(aa[i] != AminoAcidResources.getAminoAcidFromCodon(codonList[l])){
          System.err.println("AA " + aa[i] + " not found, said codon was " + codonList[l]);
        }
      }
    }
  }
  @Test
  public void AminoAcidLLTest1(){
    System.out.print("hello");
    String expected = "GAEFCHDILMNPQRVWST";
    AminoAcidLL thingy = AminoAcidLL.createFromRNASequence("GGGGCCGAGUUCUGCCACGACAUACUCAUGAACCCCCAGCGUGUGUGGAGCACGUAG");
    for (int i = 0; i < expected.length(); i++) {
      assertEquals(expected.charAt(i), thingy.aminoAcid);
      thingy =thingy.next;
    }
  }
  @Test
  public void AminoAcidLLTest2(){
    System.out.print("hello");
    AminoAcidLL thing2 = AminoAcidLL.createFromRNASequence("GGGGCCGAGUUCUGCCACGACAUACUCAUGAACCCCCAGCGUGUGUGGAGCACGUAG");
    assertEquals(false, thing2.isSorted());
  }
  @Test
  public void AminoAcidLLTest3(){
    System.out.print("hello");
    AminoAcidLL thing3 = AminoAcidLL.createFromRNASequence("ACGUAGUGG");
    assertEquals(true, thing3.isSorted());
  }
  @Test
  public void AminoAcidLLTest4(){
    System.out.print("hello");
    AminoAcidLL thing4 = AminoAcidLL.createFromRNASequence("GGGGCCGAGUUCUGCCACGACAUACUCAUGAACCCCCAGCGUGUGUGGAGCACGUAG");
    thing4 = AminoAcidLL.sort(thing4);
    assertEquals(true, thing4.isSorted());
  }
  @Test
  public void AminoAcidLLTest5(){
    System.out.print("hello");
    char[] expected = {'G','A','E','F','C','H','D','I','L','M','N','P','Q','R','V','W','S','T'};
    AminoAcidLL thing5 = AminoAcidLL.createFromRNASequence("GGGGCCGAGUUCUGCCACGACAUACUCAUGAACCCCCAGCGUGUGUGGAGCACGUAG");
    assertArrayEquals(expected, thing5.aminoAcidList());
  }
  @Test
  public void AminoAcidLLTest6(){
    System.out.print("hello");
    int[] expected = {3, 2, 1};
    String testSequence = "GGAGGAGGAGAAGAAGCUUAG";
    AminoAcidLL thing6 = AminoAcidLL.createFromRNASequence(testSequence);
    assertArrayEquals(expected, thing6.aminoAcidCounts());
  }
  @Test
  public void AminoAcidLLTest7(){
    System.out.print("hello");
    AminoAcidLL list1 = AminoAcidLL.createFromRNASequence("GAGGAGGAGACCACCUGCGACGACUAG");
    AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("GGUGGUGGUGAGGAGGAGACCACCGACUAG");
    list1 = AminoAcidLL.sort(list1);
    list2 = AminoAcidLL.sort(list2);
    assertEquals(5, list1.aminoAcidCompare(list2));
  }
  @Test
  public void AminoAcidLLTest8(){
    System.out.print("hello");
    AminoAcidLL list1 = AminoAcidLL.createFromRNASequence("GGGGAGGCGGUGUGA");
    AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("GGUGGUGGUGAGGAGGAGACCACCGACUAG");
    list1 = AminoAcidLL.sort(list1);
    list2 = AminoAcidLL.sort(list2);
    assertEquals(9, list1.aminoAcidCompare(list2));
  }
  @Test
  public void AminoAcidLLTest9(){
    System.out.print("hello");
    AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("UAG");
    AminoAcidLL list1 = AminoAcidLL.createFromRNASequence("GGUGGUGGUGAGGAGGAGACCACCGACUAG");
    list1 = AminoAcidLL.sort(list1);
    list2 = AminoAcidLL.sort(list2);
    assertEquals(9, list1.aminoAcidCompare(list2));
  }
  @Test
  public void AminoAcidLLTest10(){
    System.out.print("hello");
    AminoAcidLL list1 = AminoAcidLL.createFromRNASequence("GGGGAGGCCGCUUAG");
    AminoAcidLL list2 = AminoAcidLL.createFromRNASequence("GGAGGCGAAGCGGCAUAG");
    list1 = AminoAcidLL.sort(list1);
    list2 = AminoAcidLL.sort(list2);
    assertEquals(9, list1.codonCompare(list2));
  }
}
