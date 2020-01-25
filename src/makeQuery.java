import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

public class makeQuery{
    public static void main(String args[]){
        String sequence = "";
        double[] time = new double[100000];
        long[] memory = new long[100000];
        double[] probabilities = new double[0];
        try {
            File file = new File("chr22.maf.ancestors.42000000.complete.boreo.fa.txt");
            BufferedReader br = new BufferedReader(new FileReader(file));
            String str;
            while((str = br.readLine()) != null){
                sequence = str;
            }
            br.close();

            file = new File("chr22.maf.ancestors.42000000.complete.boreo.conf.txt");
            br = new BufferedReader((new FileReader(file)));
            String[] arr = new String[0];
            while((str = br.readLine()) != null){
                arr = str.split(" ");
            }
            probabilities = new double[arr.length];
            for(int i = 0; i < arr.length; i++){
                probabilities[i] = Double.parseDouble(arr[i]);
            }
            br.close();
            }catch (IOException e) {
            e.printStackTrace();
        }
        int correctAlignments = 0;
        Random random = new Random();
        for(int i = 0; i < 1000; i++){
            int index = random.nextInt(sequence.length()-Integer.parseInt(args[0]) - 6);
            StringBuilder query = new StringBuilder();
            int indelIndex = random.nextInt(Integer.parseInt(args[0])-6);
            int indelIndex2 =  random.nextInt(Integer.parseInt(args[0])-6);

            for(int j = index; j < index + Integer.parseInt(args[0]) + 6; j++){
                if(Math.random() > probabilities[j]){
                    query.append(randomChar(sequence.charAt(j)));
                }else{
                    query.append(sequence.charAt(j));
                }
                if(j-index == indelIndex || j-index == indelIndex2){
                    j += 6;

                }
            }
            System.out.println(index);
            System.out.println(query);
            long beforeUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
            long startTime = System.nanoTime();
            ArrayList<Alignment> alignments = findQuery.find(query.toString());
            long AfterUsedMem=Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory();
            time[i] = System.nanoTime()-startTime;
            memory[i] = beforeUsedMem - AfterUsedMem;

            boolean flag = true;
            for(Alignment align: alignments){
                if(((Math.abs(align.genomeRangeAignment[0] - index)) < 2 || (Math.abs(align.genomeRangeAignment[1] - (index + Integer.parseInt(args[0]))) < 2)) && flag) {
                    correctAlignments++;
                    System.out.println(correctAlignments);
                    System.out.println(align);
                    flag = false;
                }
            }
            System.out.println("-------------------------------");
        }
        System.out.println("Correct alignments: " + correctAlignments);
        double averageTime = 0;
        long averageMemory = 0;
        for(int i = 0; i < time.length; i++){
            averageTime += time[i];
            averageMemory += memory[i];
        }
        System.out.println("Average time: " + averageTime/1000);
        System.out.println("Average memory: " + averageMemory/1000);

    }

    public static char randomChar(char exclude){
        HashSet<Character> nucletoides = new HashSet<>();
        nucletoides.add('A');
        nucletoides.add('C');
        nucletoides.add('T');
        nucletoides.add('G');
        if(exclude != ' ') {
            nucletoides.remove(exclude);
        }
        int size = nucletoides.size();
        int item = new Random().nextInt(size); // In real life, the Random object should be rather more shared than this
        int i = 0;
        for(Character obj : nucletoides) {
            if (i == item)
                return obj;
            i++;
        }
        return 0;
    }
}
