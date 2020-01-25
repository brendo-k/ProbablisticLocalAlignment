import java.io.*;
import java.util.*;

public class bwt {
    static String sequence = "";
    public static void main(String[] args){
        long startTime = System.nanoTime();
        File file = new File(args[0]);
        String[] arr;
        double[] probabilities = null;

        try {
            BufferedReader br = new BufferedReader(new FileReader(file));
            String str;
            while((str = br.readLine()) != null){
                sequence = str;
            }
            br.close();
            arr = new String[sequence.length() + 1];
            file = new File(args[1]);
            br = new BufferedReader((new FileReader(file)));
            while((str = br.readLine()) != null){
                arr = str.split(" ");
            }
            probabilities = new double[arr.length + 1];
            for(int i = 0; i < arr.length; i++){
                probabilities[i] = Double.parseDouble(arr[i]);
            }
            br.close();
        }catch (IOException e) {
            e.printStackTrace();
        }
        sequence += '$';
        Integer[] seqStart = new Integer[sequence.length()];
        for(int i = 0; i < seqStart.length; i++){
            seqStart[i] = i;
        }

        long beforeUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
        Arrays.sort(seqStart, new Comparator(){
            public int compare(Object obj1, Object obj2){
                int index1 = (int)obj1;
                int index2 = (int)obj2;
                while(sequence.charAt((int)index1) == sequence.charAt((int)index2)){
                    index1 = (index1 + 1) % sequence.length();
                    index2 = (index2 + 1) % sequence.length();
                }
                if( sequence.charAt(index1) < sequence.charAt(index2)){
                    return -1;
                }else{
                    return 1;
                }
            }
        });
        long afterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
        System.out.println(beforeUsedMem-afterUsedMem);

        file = new File(args[2]);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(file));
            for(int i = 0; i < seqStart.length; i++){
                char line;
                double prob;
                if(seqStart[i]-1 == -1){
                    prob = 1;
                    line = '$';
                }else{
                    prob = probabilities[seqStart[i] - 1];
                    line = sequence.charAt((seqStart[i]-1));
                }
                bw.write(line + "\t" + prob);
                bw.newLine();
            }
            bw.close();
        }catch(IOException e) {
            e.printStackTrace();
        }



        int[][] oArray= new int[sequence.length()+1][5];
        for(int i = 0; i < seqStart.length; i++){
            int charIndex;
            if(seqStart[i] - 1 == -1){
                charIndex = sequence.length()-1;
            }else{
                charIndex = seqStart[i] - 1;
            }
            oArray[i+1][0] = oArray[i][0];
            oArray[i+1][1] = oArray[i][1];
            oArray[i+1][2] = oArray[i][2];
            oArray[i+1][3] = oArray[i][3];
            oArray[i+1][4] = oArray[i][4];
            switch (sequence.charAt(charIndex)) {
                case '$':
                    oArray[i+1][0] = oArray[i][0] + 1;
                    break;
                case 'A':
                    oArray[i+1][1] = oArray[i][1] + 1;
                    break;
                case 'C':
                    oArray[i+1][2] = oArray[i][2] + 1;
                    break;
                case 'G':
                    oArray[i+1][3] = oArray[i][3] + 1;
                    break;
                case 'T':
                    oArray[i+1][4] = oArray[i][4] + 1;
                    break;

            }
        }

        int[] count = new int[5];
        for(int i = 0; i < sequence.length(); i++){
            switch (sequence.charAt(seqStart[i])){
                case('$'):
                    count[0]++;
                    break;
                case('A'):
                    count[1]++;
                    break;
                case('C'):
                    count[2]++;
                    break;
                case('G'):
                    count[3]++;
                    break;
                case('T'):
                    count[4]++;
                    break;
            }
        }

        int[] cArray = new int[6];
        cArray[0] = 0;
        cArray[1] = count[0];
        cArray[2] = cArray[1] + count[1];
        cArray[3] = cArray[2] + count[2];
        cArray[4] = cArray[3] + count[3];
        cArray[5] = seqStart.length-1;


        try {
            //making occurrence array to do FM indexing and exact match
            file = new File(args[3]);
            BufferedWriter br = new BufferedWriter(new FileWriter(file));
            for (int[] ints : oArray) {
                String str = Arrays.toString(ints);
                String[] strSplit = str.substring(1, str.length() - 1).split(",");
                for(int i = 0; i < strSplit.length; i++){
                    strSplit[i] = strSplit[i].trim();
                }
                String line = String.join("\t", strSplit);
                br.write(line + "\n");
            }
            br.close();

            //making C array to index where each char.
            file = new File(args[4]);
            br = new BufferedWriter(new FileWriter(file));
            String str = Arrays.toString(cArray);
            String[] strSplit = str.substring(1, str.length() - 1).split(",");
            for(int i = 0; i < strSplit.length; i++){
                strSplit[i] = strSplit[i].trim();
            }
            String line = String.join("\t", strSplit);
            br.write(line + "\n");
            br.close();


            br = new BufferedWriter(new FileWriter(new File("start.txt")));
            for(int i = 0; i < seqStart.length; i ++){
                br.write(seqStart[i] + "\n");

            }
            br.close();

        } catch (IOException e) {
            e.printStackTrace();
        }


        System.out.println((System.nanoTime() - startTime)/1000000 );

    }
}
