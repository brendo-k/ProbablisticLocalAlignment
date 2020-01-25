import java.io.*;
import java.lang.reflect.Array;
import java.util.*;

public class findQuery {

    public static ArrayList<Alignment> find(String query){
        String bwt = "";
        String sequence = "";
        int[] cArray = null;
        int[][] oArray = null;
        Double[] probability = new Double[0];
        int[] start = new int[0];

        try{
            File file = new File("out.txt");
            BufferedReader br = new BufferedReader(new FileReader(file));
            StringBuilder sb = new StringBuilder();
            ArrayList<Double> probs = new ArrayList<>();
            String line;
            while((line = br.readLine()) != null){
                String[] str = line.split("\\t");
                probs.add(Double.parseDouble(str[1]));
                sb.append(str[0]);
            }
            bwt = sb.toString();
            br.close();
            probability =  probs.toArray(probability);

            file = new File("o.txt");
            br = new BufferedReader(new FileReader(file));
            oArray = new int[bwt.length()+1][5];
            int index = 0;
            while((line = br.readLine()) != null){
                String[] stringArr = line.split("\\t");
                for(int i = 0; i < stringArr.length; i++){
                    oArray[index][i] = Integer.parseInt(stringArr[i]);
                }
                index++;
            }
            br.close();

            file = new File("c.txt");
            br = new BufferedReader(new FileReader(file));
            cArray = new int[6];
            while((line = br.readLine()) != null){
                String[] strArr = line.split("\\t");
                for(int i = 0; i < strArr.length; i++){
                    cArray[i] = Integer.parseInt(strArr[i]);
                }
            }
            br.close();

            file = new File("start.txt");
            start = new int[bwt.length()];
            br = new BufferedReader(new FileReader(file));
            int counter = 0;
            while((line = br.readLine()) != null){
                start[counter] = Integer.parseInt(line);
                counter++;
            }
            br.close();

            file = new File("chr22.maf.ancestors.42000000.complete.boreo.fa.txt");
            br = new BufferedReader(new FileReader(file));
            sequence = br.readLine();

        }catch(IOException e){
            e.printStackTrace();
        }


        int index = query.length();
        ArrayList<State> seeds = new ArrayList<>();
        while(index > 10){
           String substring = query.substring(0, index);
            ArrayList<State> temp = match(substring, cArray, oArray, probability, bwt);
            seeds.addAll(temp);

            int maxInt = 0;
            for(State seed: temp){
                if(seed.index > maxInt){
                    maxInt = seed.index;
                }
            }
            index = maxInt;
            temp.clear();
            if(index == 0){
                break;
            }

        }

        ArrayList<Alignment> alignments = new ArrayList<>();
        for (State seed: seeds) {
            int[] rangeSequence = new int[2];
            int[] rangeQuery = new int[2];
            rangeQuery[1] = seed.startIndex;
            rangeQuery[0] = seed.index;
            for(int i = seed.lowerRange; i < seed.upperRange; i++){
                rangeSequence[1] = start[i] + (seed.startIndex - seed.index);
                rangeSequence[0] = start[i];
                Alignment align = new Alignment(sequence, rangeSequence, query, rangeQuery, probability);
                alignments.add(align);
            }
        }

//        for(Alignment align: alignments){
//            System.out.println(align);
//        }
        return(alignments);

    }

    public static ArrayList<State> match(String query, int[] cArray, int[][] oArray, Double[] probability, String bwt){
        int subsAllowed = 5;
        int minExact = 5;
        int queryIndex = query.length()-1;
        int lowPointer = 0;
        int upperPointer = cArray[5]+1;
        PriorityQueue<State> moves = new PriorityQueue<>();

        int startIndex = queryIndex;
        int numSubs = 0;
        char character = query.charAt(queryIndex);
        State previous =  new State( -1, startIndex, 0, cArray[5]+1, 'X');
        ArrayList<State> seeds = new ArrayList<>();

        while(lowPointer < upperPointer && queryIndex >= 0){
            //updating ranges
            int[] range = occurrence(oArray, lowPointer, upperPointer, character);
            lowPointer = cArray[charToInt(character)] + range[0];
            upperPointer = cArray[charToInt(character)] + range[1];
            State state = new State(queryIndex, startIndex, lowPointer, upperPointer, character);
            queryIndex--;
            if(queryIndex < 0){
                break;
            }else {
                character = query.charAt(queryIndex);
            }
            state.characterNext.add(character);
            state.getLowestProb(probability, bwt);

            if(lowPointer >= upperPointer && numSubs < subsAllowed){
                seeds.add(previous);
                numSubs++;
                if(moves.size() == 0){
                   return seeds;
                }
                state = moves.poll();
                while(state.character != 'X'){
                    if(moves.size() == 0){
                        return seeds;
                    }
                    state = moves.poll();
                }
                lowPointer = state.lowerRange;
                upperPointer = state.upperRange;
                character = state.lowestChar;
                queryIndex = state.index - 1;
                state.characterNext.add(character);
                state.getLowestProb(probability, bwt);
            }else if( query.length() - queryIndex > minExact) {
                moves.add(state);
                previous = state;
            }
        }
        if(numSubs != subsAllowed){
            seeds.add(previous);
        }

        return seeds;
    }


    public static int[] occurrence(int[][] oArray, int lowerPointer, int upperPointer, char c){
        int[] range = new int[2];
        range[0] = oArray[lowerPointer][charToInt(c)];
        range[1] = oArray[upperPointer][charToInt(c)];
        return range;
    }

    public static int charToInt(char c) {
        switch (c){
            case '$':
                return 0;
            case 'A':
                return 1;
            case 'C':
                return 2;
            case 'G':
                return 3;
            case 'T':
                return 4;
        }
        return -1;
    }



}
class State implements Comparable<State> {
    public double lowestProbability;
    public char lowestChar;
    public int lowerRange;
    public int upperRange;
    public int startIndex;
    public int index;
    public char character;
    public ArrayList<Character> characterNext;


    public State(int index, int startIndex, int lowerRange, int upperRanger, char character){
        this.index = index;
        this.startIndex = startIndex;
        this.lowerRange = lowerRange;
        this.upperRange = upperRanger;
        this.character = character;
        characterNext = new ArrayList<>();
    }

    public int compareTo(State state){
        return Double.compare(this.lowestProbability, state.lowestProbability);
    }
    public void getLowestProb(Double[] prob, String bwt){
        double lowestProb = Double.MAX_VALUE;

        for(int i = lowerRange; i < upperRange; i++){
            if(prob[i] < lowestProb && !characterNext.contains(bwt.charAt(i))){
               this.lowestProbability = prob[i];
               lowestChar = bwt.charAt(i);
            }
        }
        if(this.lowestProbability == 0){
            this.lowestProbability = 2.0;
            this.lowestChar = 'X';
        }

    }

    public String toString(){
        return Integer.toString(lowerRange) + "-" + Integer.toString(upperRange);
    }

}
