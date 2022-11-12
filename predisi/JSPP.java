import java.io.*;
import jspp.*;


public class JSPP {
  public JSPP() {
  }


  public static void main(String[] args) {
    if(args.length!=3){
      System.out.println(
          "USAGE: JSPP positive-matrix.smx sequences.fasta result.txt");
      System.exit(0);
    }
    File pmatrix=new File(args[0]);
    
    File sequences=new File(args[1]);
    File result=new File(args[2]);
    try {
      //import matrix
      BufferedReader se = new BufferedReader(new FileReader(sequences));
      ObjectInputStream in = new ObjectInputStream(new FileInputStream(pmatrix));
      SearchMatrix smp = (SearchMatrix) in.readObject();
      in.close();
       
      //construct a SignalPeptidePredictor object 
      SignalPeptidePredictor pd = new SignalPeptidePredictor(smp);
      BufferedWriter out = new BufferedWriter(new FileWriter(result));
      
      String line;
      while ( (line = se.readLine()) != null) {
        if (line.indexOf(">") == -1) {
          continue;
        }
        String id = line.substring(line.indexOf(">") + 1);
        String sequence;
        if ( (sequence = se.readLine()) == null) {
          System.err.println("Error reading FASTA-File !");
          System.exit(0);
        }

        int cpos = pd.predictEnhancedPosition(sequence);
        String yn;
        if (pd.isSignalPeptide()) {
          yn = "Y";
        }
        else {
          yn = "N";
        }
        double score=pd.getScore();
        out.write(id + "\t" + cpos + "\t" + yn + "\t"+score+"\n");
      }
      se.close();
      out.close();
    }
    catch (IOException ex) {
      System.err.println("IO-Error: "+ex);
    }
    catch (ClassNotFoundException ex) {
      System.err.println("Error reading matrices (wrong file ?)");
    }


  }

}