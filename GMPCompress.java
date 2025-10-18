package gmp;

import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import org.tukaani.xz.LZMA2Options;
import org.tukaani.xz.XZOutputStream;

public class GMPCompress {
    // Maximum code table size (PrCT size)
    private static final int maxPrCTSize = 4;
    // Population size for genetic algorithm (will be set dynamically based on unique characters)
    private static int Pop;
    // Maximum generations per evolution cycle
    private static final int maxGenerations = 10;
    
    public static void main(String[] args) {
        // Set your dataset name here (e.g., "HT")
        String datasetName = "PDBaa";
        Path cwd = Paths.get("").toAbsolutePath();
        Path datasetPath = cwd.resolve(datasetName);
        Path patternsOutPath = cwd.resolve(datasetName + "PrCT.txt");
        Path dictionaryOutPath = cwd.resolve(datasetName + "dictionary.xz");
        Path updatedOutPath = cwd.resolve(datasetName + "update");

        try {
            if (!Files.exists(datasetPath)) {
                System.err.println("Dataset file not found: " + datasetPath);
                System.exit(2);
            }

            String dataset = Files.readString(datasetPath, StandardCharsets.UTF_8);
            dataset = dataset.replace("\r\n", "\n");

            // 1) Collect unique characters present in dataset (AA) - including non-standard amino acids
            List<Character> aminoAcids = collectUniqueCharacters(dataset);
            if (aminoAcids.isEmpty()) {
                System.err.println("Dataset appears empty or contains no characters.");
                System.exit(3);
            }
            
            // Set Pop to the number of unique characters found
            Pop = aminoAcids.size();
            System.out.println("Found " + Pop + " unique characters. Setting population size to " + Pop);
            
            // Save unique characters to file
            Path uniqueCharsPath = cwd.resolve(datasetName + "_unique_chars.txt");
            saveUniqueCharacters(aminoAcids, uniqueCharsPath);

            // 2) Initialize population Pop with random kAAMs from amino acid alphabet AA
            Random rng = new Random();
            List<String> population = initializePopulation(aminoAcids, rng);
            
            // 3) Initialize PrCT (Pattern Code Table)
            List<String> PrCT = new ArrayList<>();

            // 4) Main GMP algorithm loop: while |PrCT| < maxPrCTSize
            while (PrCT.size() < maxPrCTSize) {
                System.out.println("Building pattern " + (PrCT.size() + 1) + "/" + maxPrCTSize);
                
                // Select P1, P2 from Pop (biased by usage/compression gain)
                String[] parents = selectParents(population, dataset, PrCT, rng);
                String P1 = parents[0];
                String P2 = parents[1];

                // Evolution loop: repeat until max generations
                for (int gen = 0; gen < maxGenerations; gen++) {
                    // Crossover: C1, C2 ← Crossover(P1, P2)
                    String[] children = singlePointCrossover(P1, P2, rng);
                    String C1 = children[0];
                    String C2 = children[1];

                    // Mutation: C1 ← Mutation(C1, Pop), C2 ← Mutation(C2, Pop)
                    C1 = mutate(C1, aminoAcids, rng);
                    C2 = mutate(C2, aminoAcids, rng);

                    // Evaluate fitness and update parents if children are better
                    if (fitness(C1, dataset, PrCT) > fitness(P1, dataset, PrCT)) {
                        P1 = C1;
                    }
                    if (fitness(C2, dataset, PrCT) > fitness(P2, dataset, PrCT)) {
                        P2 = C2;
                    }
                }

                // Add best candidates {P1, P2} to PrCT
                PrCT.add(P1);
                if (PrCT.size() < maxPrCTSize) {
                    PrCT.add(P2);
                }
                
                // Update population with new patterns
                updatePopulation(population, P1, P2, rng);
            }

            // Save evolved patterns (PrCT)
            try (BufferedWriter w = Files.newBufferedWriter(patternsOutPath, StandardCharsets.UTF_8)) {
                for (String pattern : PrCT) {
                    w.write(pattern);
                    w.newLine();
                }
            }

            // 5) Build dictionary D: pattern → positions from PrCT
            String working = dataset;

            try (OutputStream fileOut = Files.newOutputStream(dictionaryOutPath);
                 BufferedOutputStream bufferedOut = new BufferedOutputStream(fileOut)) {
                LZMA2Options options = new LZMA2Options();
                options.setPreset(LZMA2Options.PRESET_MAX);
                try (XZOutputStream xzOut = new XZOutputStream(bufferedOut, options)) {
                    for (String pattern : PrCT) {
                        List<Integer> positions = new ArrayList<>();
                        int fromIndex = 0;
                        int patternLength = pattern.length();
                        String dots = repeatDot(patternLength);
                        while (fromIndex <= working.length() - patternLength) {
                            int idx = working.indexOf(pattern, fromIndex);
                            if (idx == -1) break;
                            positions.add(idx);
                            StringBuilder sb = new StringBuilder(working);
                            sb.replace(idx, idx + patternLength, dots);
                            working = sb.toString();
                            fromIndex = idx + 1; // allow overlap
                        }

                        StringBuilder line = new StringBuilder();
                        line.append(pattern);
                        if (!positions.isEmpty()) {
                            int prev = 0;
                            for (int i = 0; i < positions.size(); i++) {
                                int absolute = positions.get(i);
                                int delta = (i == 0) ? absolute : (absolute - prev);
                                line.append(' ').append(delta);
                                prev = absolute;
                            }
                        }
                        line.append('\n');
                        xzOut.write(line.toString().getBytes(StandardCharsets.UTF_8));
                    }
                    xzOut.finish();
                }
            }

            // 6) Filter PD by excluding recorded patterns, yielding D'
            String updated = working.replace(".", "");
            try (BufferedWriter writer = Files.newBufferedWriter(updatedOutPath, StandardCharsets.UTF_8)) {
                writer.write(updated);
            }

            System.out.println("GMP Compression completed!");
            System.out.println("Population size (Pop): " + Pop);
            System.out.println("Saved unique characters (" + aminoAcids.size() + " chars): " + uniqueCharsPath.toAbsolutePath());
            System.out.println("Saved PrCT (" + PrCT.size() + " patterns): " + patternsOutPath.toAbsolutePath());
            System.out.println("Saved compressed dictionary (XZ): " + dictionaryOutPath.toAbsolutePath());
            System.out.println("Saved updated database: " + updatedOutPath.toAbsolutePath());
        } catch (IOException e) {
            System.err.println("I/O error: " + e.getMessage());
            System.exit(4);
        }
    }

    // Initialize population Pop with random kAAMs from amino acid alphabet AA
    private static List<String> initializePopulation(List<Character> aminoAcids, Random rng) {
        List<String> population = new ArrayList<>();
        for (int i = 0; i < Pop; i++) {
            population.add(randomPattern(aminoAcids, rng));
        }
        return population;
    }

    // Select P1, P2 from Pop (biased by usage/compression gain)
    private static String[] selectParents(List<String> population, String dataset, List<String> PrCT, Random rng) {
        // Tournament selection with bias towards better fitness
        String P1 = tournamentSelection(population, dataset, PrCT, rng);
        String P2 = tournamentSelection(population, dataset, PrCT, rng);
        
        // Ensure P1 and P2 are different
        while (P2.equals(P1) && population.size() > 1) {
            P2 = tournamentSelection(population, dataset, PrCT, rng);
        }
        
        return new String[]{P1, P2};
    }

    // Tournament selection for parent selection
    private static String tournamentSelection(List<String> population, String dataset, List<String> PrCT, Random rng) {
        int tournamentSize = Math.min(3, population.size());
        String best = population.get(rng.nextInt(population.size()));
        double bestFitness = fitness(best, dataset, PrCT);
        
        for (int i = 1; i < tournamentSize; i++) {
            String candidate = population.get(rng.nextInt(population.size()));
            double candidateFitness = fitness(candidate, dataset, PrCT);
            if (candidateFitness > bestFitness) {
                best = candidate;
                bestFitness = candidateFitness;
            }
        }
        
        return best;
    }

    // Update population with new patterns
    private static void updatePopulation(List<String> population, String P1, String P2, Random rng) {
        // Replace worst individuals with new patterns
        population.remove(rng.nextInt(population.size()));
        population.add(P1);
        
        if (population.size() < Pop) {
            population.remove(rng.nextInt(population.size()));
            population.add(P2);
        }
    }

    // Fitness function based on compression improvement
    private static double fitness(String pattern, String dataset, List<String> PrCT) {
        // Count occurrences of the pattern
        long occurrences = scorePattern(dataset, pattern);
        
        if (occurrences == 0) return 0;
        
        // Fitness is proportional to compression gain
        // More occurrences and longer patterns give better compression
        double compressionGain = occurrences * pattern.length();
        
        // Penalize patterns that are too similar to existing ones
        for (String existingPattern : PrCT) {
            if (pattern.contains(existingPattern) || existingPattern.contains(pattern)) {
                compressionGain *= 0.5; // Reduce fitness for similar patterns
            }
        }
        
        return compressionGain;
    }

    private static List<Character> collectUniqueCharacters(String text) {
        Set<Character> set = new HashSet<>();
        for (int i = 0; i < text.length(); i++) {
            char c = text.charAt(i);
            if (c == '\n' || c == '\r') continue;
            set.add(c);
        }
        return new ArrayList<>(set);
    }

    private static void saveUniqueCharacters(List<Character> characters, Path filePath) throws IOException {
        try (BufferedWriter writer = Files.newBufferedWriter(filePath, StandardCharsets.UTF_8)) {
            for (Character c : characters) {
                writer.write(c);
                writer.newLine();
            }
        }
    }

    private static String randomPattern(List<Character> aminoAcids, Random rng) {
        int len = 2 + rng.nextInt(3); // 2..4
        StringBuilder sb = new StringBuilder(len);
        for (int i = 0; i < len; i++) {
            char c = aminoAcids.get(rng.nextInt(aminoAcids.size()));
            sb.append(c);
        }
        return sb.toString();
    }

    private static String[] singlePointCrossover(String p1, String p2, Random rng) {
        int s = Math.min(p1.length(), p2.length());
        int cp = 1 + rng.nextInt(s); // 1..s
        String c1 = p1.substring(0, cp) + p2.substring(Math.min(cp, p2.length()));
        String c2 = p2.substring(0, cp) + p1.substring(Math.min(cp, p1.length()));
        return new String[] { c1, c2 };
    }

    private static String mutate(String pattern, List<Character> aminoAcids, Random rng) {
        int i = rng.nextInt(pattern.length());
        char newAA = aminoAcids.get(rng.nextInt(aminoAcids.size()));
        StringBuilder sb = new StringBuilder(pattern);
        sb.setCharAt(i, newAA);
        return sb.toString();
    }

    private static long scorePattern(String text, String pattern) {
        long count = 0;
        int from = 0;
        int len = pattern.length();
        while (from <= text.length() - len) {
            int idx = text.indexOf(pattern, from);
            if (idx == -1) break;
            count++;
            from = idx + 1; // allow overlap
        }
        return count;
    }

    private static String repeatDot(int count) {
        if (count <= 0) return "";
        StringBuilder sb = new StringBuilder(count);
        for (int i = 0; i < count; i++) sb.append('.');
        return sb.toString();
    }
}



