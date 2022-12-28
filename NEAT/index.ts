const maxNodes = 1000000;


class RandomHashSet<TYPE> {
    constructor() {
        this.hashSet = new Set();
        this.arrayList = [];
    }

    hashSet: Set<TYPE>;
    arrayList: TYPE[];

    contains(object: TYPE): boolean {
        return this.hashSet.has(object);
    }

    randomElement(): TYPE | null {
        return this.arrayList.length > 0 ? this.arrayList[Math.floor(Math.random() * this.arrayList.length)] : null;
    }

    size(): number {
        return this.arrayList.length;
    }

    add(object: TYPE): void {
        if (!this.hashSet.has(object)) {
            this.hashSet.add(object);
            this.arrayList.push(object);
        }
    }

    addSorted(object: any) {
        for (let i = 0; i < this.size(); i++) {
            const innovation = (this.arrayList[i] as Gene).getInnovationNumber();

            if (object.getInnovationNumber() < innovation) {
                this.arrayList[i] = object;
                this.hashSet.add(object);
                return;
            }

        }

        this.add(object)
    }

    clear(): void {
        this.hashSet.clear();
        this.arrayList = [];
    }

    get(index: number): TYPE | null {
        if (index < 0 || index >= this.arrayList.length) return null;
        return this.arrayList[index];
    }

    remove(index: number): void {
        if (index < 0 || index >= this.arrayList.length) return;
        this.hashSet.delete(this.arrayList[index]);
        this.arrayList.splice(index, 1);
    }

    removeElement(object: TYPE): void {
        if (this.hashSet.has(object)) {
            this.hashSet.delete(object);
            this.arrayList.splice(this.arrayList.indexOf(object), 1);
        }
    }

    getData(): TYPE[] {
        return this.arrayList;
    }
}

class RandomSelector<TYPE> {
    constructor() {
        this.objects = [];
        this.scores = [];
    }

    objects: TYPE[] = [];
    scores: number[] = [];
    total_score = 0;

    add(object: TYPE, score: number): void {
        this.objects.push(object);
        this.scores.push(score);
        this.total_score += score;
    }

    random(): TYPE | null {
        if (this.objects.length == 0) return null;

        let random = Math.random() * this.total_score;
        let c = 0;
        for (let i = 0; i < this.objects.length; i++) {
            c += this.scores[i];
            if (c > random) return this.objects[i];
        }

        return null;
    }

    public reset(): void {
        this.objects = [];
        this.scores = [];
        this.total_score = 0;
    }
}

class Gene {
    constructor(innovationNumber: number) {
        this.innovationNumber = innovationNumber;
    }

    innovationNumber: number;

    getInnovationNumber(): number {
        return this.innovationNumber;
    }

    setInnovationNumber(innovationNumber: number): void {
        this.innovationNumber = innovationNumber;
    }
}

class NodeGene extends Gene {
    constructor(innovationNumber: number) {
        super(innovationNumber);
        this.x = 0;
        this.y = 0;
    }

    type: "" | "input" | "output" = "";

    getType(): "" | "input" | "output" {
        return this.type;
    }

    setType(type: "" | "input" | "output"): void {
        this.type = type;
    }

    getX(): number {
        return this.x;
    }

    getY(): number {
        return this.y;
    }

    setX(x: number): void {
        this.x = x;
    }

    setY(y: number): void {
        this.y = y;
    }

    equals(nodeGene: NodeGene): boolean {
        return this.innovationNumber == nodeGene.innovationNumber;
    }

    hashCode(): number {
        return this.innovationNumber;
    }

    x: number;
    y: number;
}

class ConnectionGene extends Gene {
    constructor(from: NodeGene, to: NodeGene) {
        super(0);
        this.from = from;
        this.to = to;
    }

    from: NodeGene;
    to: NodeGene;
    weight: number = 0;
    enabled: boolean = true;

    replaceIndex: number | undefined;

    setReplaceIndex(replaceIndex: number): void {
        this.replaceIndex = replaceIndex;
    }

    getReplaceIndex(): number | undefined {
        return this.replaceIndex;
    }

    getFrom(): NodeGene {
        return this.from;
    }

    getTo(): NodeGene {
        return this.to;
    }

    setFrom(from: NodeGene): void {
        this.from = from;
    }

    setTo(to: NodeGene): void {
        this.to = to;
    }

    getWeight(): number {
        return this.weight;
    }

    setWeight(weight: number): void {
        this.weight = weight;
    }

    isEnabled(): boolean {
        return this.enabled;
    }

    setEnabled(enabled: boolean): void {
        this.enabled = enabled;
    }

    equals(connectionGene: ConnectionGene): boolean {
        return this.from.equals(connectionGene.from) && this.to.equals(connectionGene.to);
    }

    hashCode(): number {
        return this.from.getInnovationNumber() * maxNodes + this.to.getInnovationNumber();
    }
}

class Genome {
    connections = new RandomHashSet<ConnectionGene>();
    nodes = new RandomHashSet<NodeGene>();
    neat: Neat;

    constructor(neat: Neat) {
        this.neat = neat;
    }

    distance(g2: Genome): number {

        let g1: Genome = this

        let highestInnovationG1 = (g1.getConnections().size() != 0 ? g1.getConnections().get(g1.getConnections().size() - 1)?.getInnovationNumber() : 0) ?? 0
        let highestInnovationG2 = (g2.getConnections().size() != 0 ? g2.getConnections().get(g2.getConnections().size() - 1)?.getInnovationNumber() : 0) ?? 0

        if (highestInnovationG1 < highestInnovationG2) {
            const temp = g1;
            g1 = g2;
            g2 = temp;
        }

        let indexg1 = 0;
        let indexg2 = 0;

        let disjoint = 0;
        let excess = 0;
        let weightDiff = 0;
        let similar = 0;

        while (indexg1 < g1.getConnections().size() && indexg2 < g2.getConnections().size()) {
            const c1: ConnectionGene | null = g1.getConnections().get(indexg1)
            const c2: ConnectionGene | null = g2.getConnections().get(indexg2)

            const in1: number | null = c1?.getInnovationNumber() ?? null;
            const in2: number | null = c2?.getInnovationNumber() ?? null;

            if (in1 == in2) {
                similar++;
                weightDiff += Math.abs((c1?.getWeight() ?? 0) - (c2?.getWeight() ?? 0));
                indexg1++;
                indexg2++;
            } else if ((in1 ?? 0) > (in2 ?? 0)) {
                disjoint++;
                indexg2++;
            } else {
                disjoint++;
                indexg1++;
            }
        }

        weightDiff /= similar;
        excess = g1.getConnections().size() - indexg1;

        let N = Math.max(g1.getConnections().size(), g2.getConnections().size());
        if (N < 20) N = 1;

        return neat.getC1() * disjoint / N + neat.getC2() * excess / N + neat.getC3() * weightDiff
    }

    crossover(g1: Genome, g2: Genome): Genome {

        const neat = g1.getNeat();
        const child = neat.emptyGenome();

        let indexg1 = 0;
        let indexg2 = 0;

        while (indexg1 < g1.getConnections().size() && indexg2 < g2.getConnections().size()) {
            const gene1: ConnectionGene | null = g1.getConnections().get(indexg1)
            const gene2: ConnectionGene | null = g2.getConnections().get(indexg2)

            const in1: number = (gene1?.getInnovationNumber()) ?? 0;
            const in2: number = (gene2?.getInnovationNumber()) ?? 0;

            if (in1 == in2 && gene1 && gene2) {

                if (Math.random() < 0.5) {
                    child.getConnections().add(neat.getConnection(gene1));
                } else {
                    child.getConnections().add(neat.getConnection(gene2));
                }

                indexg1++;
                indexg2++;
            } else if (in1 > in2) {
                indexg2++;
            } else {
                if (gene1) child.getConnections().add(neat.getConnection(gene1));
                indexg1++;
            }
        }

        while (indexg1 < g1.getConnections().size()) {
            const gene1: ConnectionGene | null = g1.getConnections().get(indexg1)
            if (gene1) child.getConnections().add(neat.getConnection(gene1));
            indexg1++;
        }

        child.getConnections().arrayList.forEach(connection => {
            child.getNodes().add(connection.getFrom());
            child.getNodes().add(connection.getTo());
        });

        return child;
    }

    mutate(): void {
        if (neat.getPROBABILITY_MutateNode() > Math.random()) this.mutateNode();
        if (neat.getPROBABILITY_MutateLink() > Math.random()) this.mutateLink();
        if (neat.getPROBABILITY_MutateToggleLink() > Math.random()) this.mutateLinkToggle();
        if (neat.getPROBABILITY_MutateWeightRandom() > Math.random()) this.mutateWeightRandom();
        if (neat.getPROBABILITY_MutateWeightShift() > Math.random()) this.mutateWeightShift();
    }

    mutateLink(): void {
        for (let i = 0; i < 100; i++) {
            const a = this.nodes.randomElement();
            const b = this.nodes.randomElement();

            if (a?.getX() == b?.getX()) continue;

            let connection
            if ((a?.getX() ?? 0) < (b?.getX() ?? 0)) {
                connection = new ConnectionGene(a!, b!);
            } else {
                connection = new ConnectionGene(b!, a!);
            }

            if (this.connections.contains(connection)) continue;

            connection = neat.getConnection(connection.getFrom(), connection.getTo());
            connection.setWeight((Math.random() * 2 - 1) * neat.getWeightRandomStrength());

            this.connections.addSorted(connection);
            return;
        }
    }

    mutateNode(): void {
        const connection = this.connections.randomElement();
        if (!connection?.getFrom) return;

        const from = connection.getFrom()
        const to = connection.getFrom()

        const replaceIndex = neat.getReplaceIndex(from, to)

        let middle;

        if (replaceIndex == 0) {
            let middle = neat.getNode()
            middle?.setX((from.getX() + to.getX()) / 2)
            middle?.setY((from.getY() + to.getY()) / 2 * Math.random() * 0.1 - 0.05)
            neat.setReplaceIndex(from, to, middle?.getInnovationNumber() ?? 0)
        } else {
            middle = neat.getNode(replaceIndex)
        }

        if (!middle) return;

        const con1 = neat.getConnection(from, middle!);
        const con2 = neat.getConnection(middle!, to);

        con1.setWeight(1);
        con2.setWeight(connection.getWeight());
        con2.setEnabled(connection.isEnabled());

        this.connections.remove(this.connections.getData().indexOf(connection));
        this.connections.add(con1);
        this.connections.add(con2);

        this.nodes.add(middle!);
    }

    mutateWeightShift(): void {
        const con = this.connections.randomElement();
        if (con?.getWeight) {
            con.setWeight(con.getWeight() + (Math.random() * 2 - 1) * neat.getWeightRandomStrength());
        }
    }

    mutateWeightRandom(): void {
        const con = this.connections.randomElement();
        if (con?.setWeight) {
            con.setWeight((Math.random() * 2 - 1) * neat.getWeightRandomStrength());
        }
    }

    mutateLinkToggle(): void {
        const con = this.connections.randomElement();
        if (con?.isEnabled) {
            con.setEnabled(!con.isEnabled());
        }
    }

    getConnections(): RandomHashSet<ConnectionGene> {
        return this.connections;
    }

    getNodes(): RandomHashSet<NodeGene> {
        return this.nodes;
    }

    getNeat(): Neat {
        return this.neat;
    }
}

class Neat {

    allConnections = new Map<ConnectionGene, ConnectionGene>();
    allNodes = new RandomHashSet<NodeGene>();

    C1 = 1;
    C2 = 1;
    C3 = 1;
    CP = 4;
    survivors = 0.8;
    weightShiftStrength = 0.3;
    weightRandomStrength = 1;

    PROBABILITY_MutateLink = 0.01;
    PROBABILITY_MutateNode = 0.1;
    PROBABILITY_MutateWeightShift = 0.02;
    PROBABILITY_MutateWeightRandom = 0.02;
    PROBABILITY_MutateToggleLink = 0;

    clients = new RandomHashSet<Client>();
    species = new RandomHashSet<Species>();

    inputSize: number = 0;
    outputSize: number = 0;
    maxClients: number = 0;

    constructor(inputSize: number, outputSize: number, clients: number) {
        this.reset(inputSize, outputSize, clients);
    }

    setReplaceIndex(node1: NodeGene, node2: NodeGene, index: number) {
        this.allConnections.get(new ConnectionGene(node1, node2))?.setReplaceIndex(index);
    }

    getReplaceIndex(node1: NodeGene, node2: NodeGene): number {
        const con = new ConnectionGene(node1, node2);
        const data = this.allConnections.get(con);
        if (!data) return 0;
        return data.getReplaceIndex() ?? 0;
    }

    emptyGenome(): Genome {
        const genome = new Genome(this);

        for (let i = 0; i < this.inputSize + this.outputSize; i++) {
            genome.getNodes().add(this.getNode(i + 1)!);
        }

        return genome;
    }

    reset(inputSize: number, outputSize: number, clients: number): void {
        this.inputSize = inputSize;
        this.outputSize = outputSize;
        this.maxClients = clients;

        this.allConnections.clear();
        this.allNodes.clear();
        this.clients.clear();

        for (let i = 0; i < inputSize; i++) {
            const node = this.getNode();
            node?.setX(0.1);
            node?.setType("input");
            node?.setY((i + 1) / (inputSize + 1))
        }

        for (let i = 0; i < outputSize; i++) {
            const node = this.getNode();
            node?.setX(0.9);
            node?.setType("output");
            node?.setY((i + 1) / (outputSize + 1))
        }

        for (let i = 0; i < this.maxClients; i++) {
            const client = new Client();
            client.setGenome(this.emptyGenome());
            client.generateCalculator();
            this.clients.add(client);
        }
    }

    getClient(index: number): Client | null {
        if (index < this.clients.size()) return null;
        return this.clients.get(index);
    }

    copyConnection(connection: ConnectionGene) {
        const newConnection = new ConnectionGene(connection.getFrom(), connection.getTo());
        newConnection.setWeight(connection.getWeight());
        newConnection.setEnabled(connection.isEnabled());
        return newConnection;
    }

    getConnection(from: NodeGene | ConnectionGene, to?: NodeGene): ConnectionGene {
        if (to && from instanceof NodeGene && !(from instanceof ConnectionGene)) {
            const connectionGene = new ConnectionGene(from, to);

            const geneInConnections = this.allConnections.get(connectionGene)
            if (geneInConnections) {
                connectionGene.setInnovationNumber(geneInConnections.getInnovationNumber());
            } else {
                connectionGene.setInnovationNumber(this.allConnections.size + 1);
                this.allConnections.set(connectionGene, connectionGene);
            }

            return connectionGene
        } else {
            if (!(from instanceof ConnectionGene)) return {} as ConnectionGene; // impossible but TS won't shut up
            const connection = new ConnectionGene(from.getFrom(), from.getTo());
            connection.setInnovationNumber(from.getInnovationNumber());
            connection.setWeight(from.getWeight());
            connection.setEnabled(from.isEnabled());
            return connection;
        }
    }

    getNode(id?: number): NodeGene | null {
        if (!id) {
            const node = new NodeGene(this.allNodes.size() + 1);
            this.allNodes.add(node);
            return node;
        } else {
            if (id <= this.allNodes.size()) return this.allNodes.get(id - 1);
            else {
                return this.getNode();
            }
        }
    }

    evolve(): void {
        this.generateSpecies();
        this.kill();
        this.removeExtinctSpecies();
        this.reproduce();
        this.mutate();

        this.clients.getData().forEach(client => {
            client.generateCalculator();
        });
    }

    generateSpecies(): void {
        this.species.getData().forEach(species => {
            species.reset()
        });

        for (let i = 0; i < this.clients.getData().length; i++) {
            const client = this.clients.getData()[i];
            if (client.getSpecies() != null) continue;

            let found = false;

            for (let i = 0; i < this.species.getData().length; i++) {
                const species = this.species.getData()[i];
                if (species.put(client)) {
                    found = true;
                    break;
                };
            };

            if (!found) {
                this.species.add(new Species(client))
            }
        };

        this.species.getData().forEach(species => {
            species.evaluateScore();
        });
    }

    kill() {
        this.species.getData().forEach(species => {
            species.kill(1 - this.survivors)
        });
    }

    removeExtinctSpecies() {
        for (let i = this.species.size() - 1; i >= 0; i--) {
            if (this.species.size() <= 1) {
                this.species.get(i)?.goExtinct()
                this.species.remove(i);
            }
        }
    }

    reproduce() {
        const selector = new RandomSelector<Species>();
        this.species.getData().forEach(species => {
            selector.add(species, species.getScore());
        });

        this.clients.getData().forEach(client => {
            if (client.getSpecies() == null) {
                const species = selector.random();
                const genome = species?.breed();
                if (!genome) return;
                client.setGenome(genome);
                species?.forcePut(client);
            }
        });
    }

    mutate() {
        this.clients.getData().forEach(client => {
            client.mutate();
        });
    }

    getC1(): number {
        return this.C1;
    }

    getC2(): number {
        return this.C2;
    }

    getC3(): number {
        return this.C3;
    }

    getWeightShiftStrength(): number {
        return this.weightShiftStrength;
    }

    getWeightRandomStrength(): number {
        return this.weightRandomStrength;
    }

    getPROBABILITY_MutateLink(): number {
        return this.PROBABILITY_MutateLink;
    }

    getPROBABILITY_MutateNode(): number {
        return this.PROBABILITY_MutateNode;
    }

    getPROBABILITY_MutateWeightShift(): number {
        return this.PROBABILITY_MutateWeightShift;
    }

    getPROBABILITY_MutateWeightRandom(): number {
        return this.PROBABILITY_MutateWeightRandom;
    }

    getPROBABILITY_MutateToggleLink(): number {
        return this.PROBABILITY_MutateToggleLink;
    }

    getCP(): number {
        return this.CP;
    }

    setCP(CP: number): void {
        this.CP = CP;
    }

    printSpecies() {
        // console.log("######################################")
        // console.log("Species: ", this.species.getData())
    }
}

class Connection {
    constructor(from: _Node, to: _Node) {
        this.from = from;
        this.to = to;
    }

    from: _Node;
    to: _Node;
    weight: number = 0;
    enabled: boolean = true;

    getFrom(): _Node {
        return this.from;
    }

    getTo(): _Node {
        return this.to;
    }

    setFrom(from: _Node): void {
        this.from = from;
    }

    setTo(to: _Node): void {
        this.to = to;
    }

    getWeight(): number {
        return this.weight;
    }

    setWeight(weight: number): void {
        this.weight = weight;
    }

    isEnabled(): boolean {
        return this.enabled;
    }

    setEnabled(enabled: boolean): void {
        this.enabled = enabled;
    }
}

class _Node {
    x: number;
    output: number = 0;
    connections: Connection[] = [];

    constructor(x: number) {
        this.x = x;
    }

    calculate(): void {
        let sum = 0;
        this.connections.forEach(connection => {
            if (connection.isEnabled()) {
                sum += connection.getWeight() * connection.getFrom().getOutput();
            }
        });

        this.output = this.activationFunction(sum);
    }

    activationFunction(x: number): number {
        return 1 / (1 + Math.exp(-x));
    }

    getX(): number {
        return this.x;
    }

    setX(x: number): void {
        this.x = x;
    }

    getOutput(): number {
        return this.output;
    }

    setOutput(output: number): void {
        this.output = output;
    }

    getConnections(): Connection[] {
        return this.connections;
    }

    setConnections(connections: Connection[]): void {
        this.connections = connections;
    }

    compareTo(node: _Node): number {
        if (this.x > node.getX()) return -1;
        if (this.x < node.getX()) return 1;
        return 0;
    }
}

class Calculator {

    inputNodes: _Node[] = [];
    outputNodes: _Node[] = [];
    hiddenNodes: _Node[] = [];

    constructor(genome: Genome) {
        const nodes = genome.getNodes();
        const connections = genome.getConnections();

        const nodeHashMap = new Map();

        nodes.arrayList.forEach((n: NodeGene) => {
            const node = new _Node(n.getX())
            nodeHashMap.set(n.getInnovationNumber(), node);

            if (n.getType() === "input") {
                this.inputNodes.push(node);
            } else if (n.getType() === "output") {
                this.outputNodes.push(node);
            } else {
                this.hiddenNodes.push(node);
            }
        });

        this.hiddenNodes.sort((a, b) => a.compareTo(b));

        connections.arrayList.forEach((connection: ConnectionGene) => {
            const from = connection.getFrom();
            const to = connection.getTo();

            const nodeFrom: _Node = nodeHashMap.get(from.getInnovationNumber());
            const nodeTo: _Node = nodeHashMap.get(to.getInnovationNumber());

            const newConnection = new Connection(nodeFrom, nodeTo);
            newConnection.setWeight(connection.getWeight());
            newConnection.setEnabled(connection.isEnabled());

            nodeTo.connections.push(newConnection);
        });
    }

    calculate(input: number[]): number[] {
        for (let i = 0; i < this.inputNodes.length; i++) {
            this.inputNodes[i].setOutput(input[i]);
        }

        this.hiddenNodes.forEach(node => {
            node.calculate();
        });

        const output: number[] = new Array(this.outputNodes.length);

        for (let i = 0; i < this.outputNodes.length; i++) {
            this.outputNodes[i].calculate();
            output[i] = this.outputNodes[i].getOutput();
        }

        return output;
    }
}

class Client {

    calculator: Calculator | undefined;
    genome: Genome | null = null;
    score: number = 0;
    species: Species | null = null;

    distance(other: Client): number {
        return this.getGenome()!.distance(other.genome!);
    }

    mutate(): void {
        this.getGenome()!.mutate();
    }

    generateCalculator(): void {
        this.calculator = new Calculator(this.genome!);
    }

    calculate(input: number[]): number[] {
        if (input.length != neat.inputSize) throw new Error("Input size is not correct");
        if (!this.calculator) this.generateCalculator();
        if (this.calculator != null) {
            return this.calculator.calculate(input);
        }
        return [-1]
    }

    getCalculator(): Calculator | undefined {
        return this.calculator;
    }

    getGenome(): Genome | null {
        return this.genome;
    }

    setGenome(genome: Genome): void {
        this.genome = genome;
    }

    getScore(): number {
        return this.score;
    }

    setScore(score: number): void {
        this.score = score;
    }

    getSpecies(): Species | null {
        return this.species;
    }

    setSpecies(Species: Species | null): void {
        this.species = Species;
    }
}

class Species {
    clients = new RandomHashSet<Client>();
    representative: Client | null;
    score: number = 0;

    constructor(representative: Client) {
        this.representative = representative;
        this.representative.setSpecies(this);
        this.clients.add(representative)
    }

    put(client: Client): boolean {
        if (client.distance(this.representative!) > (this.representative?.getGenome()?.getNeat().getCP() ?? 0)) {
            client.setSpecies(this);
            this.clients.add(client);
            return true;
        };
        return false;
    }

    forcePut(client: Client) {
        client.setSpecies(this);
        this.clients.add(client);
    }

    goExtinct(): void {
        this.clients.getData().forEach(client => {
            client.setSpecies(null);
        });
    }

    evaluateScore(): void {
        let v = 0;
        this.clients.getData().forEach(client => {
            v += client.getScore();
        });

        this.score = v / this.clients.size();
    }

    reset() {
        this.representative = this.clients.randomElement();

        this.clients.getData().forEach(client => {
            client.setSpecies(null);
        });
        this.clients.clear();

        this.clients.add(this.representative!);
        this.representative?.setSpecies(this);

        this.score = 0;
    }

    kill(percentage: number): void {
        this.clients.getData().sort((a, b) => {
            if (a.getScore() > b.getScore()) return -1;
            if (a.getScore() < b.getScore()) return 1;
            return 0;
        });

        const amount = Math.round(percentage * this.clients.size());
        for (let i = 0; i < amount; i++) {
            this.clients.get(0)?.setSpecies(null);
            this.clients.remove(0);
        }
    }

    breed(): Genome | null {
        const c1 = this.clients.randomElement();
        const c2 = this.clients.randomElement();

        if (!c1?.getGenome() && !c2?.getGenome()) return null;
        else {
            if (c1!.getScore() > c2!.getScore()) return genome.crossover(c1?.getGenome(), c2?.getGenome());
            return genome.crossover(c2.getGenome(), c1.getGenome());
        }
    }

    size() {
        return this.clients.size();
    }

    getClients(): RandomHashSet<Client> {
        return this.clients;
    }

    getRepresentative(): Client | null {
        return this.representative;
    }

    getScore(): number {
        return this.score;
    }
};

let neat: Neat;
let genome: Genome;
let input: number[];
neat = new Neat(2, 1, 1000);

for (let i = 0; i < 100; i++) {
    neat.clients.getData().forEach(client => {
        input = [Math.round(Math.random()), Math.round(Math.random())]
        const output = Math.round(client.calculate(input)[0])
        const correctOutput = input[0] ^ input[1]

        const score = 1 / (1 + Math.abs(output - correctOutput));

        client.setScore(score);
    });

    neat.evolve();
    neat.printSpecies()
}

const results = new Set<Number>();

neat.clients.getData().forEach(client => {
    const output = client.calculate([1, 0])[0]

    const correctOutput = 1

    const score = 1 / (1 + Math.abs(output - correctOutput));
    results.add(Math.round(score))
    client.setScore(score);
});

console.log(results)