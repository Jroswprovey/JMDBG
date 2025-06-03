package org.conncoll;


import java.util.function.Consumer;

public class Command {
    private final String name;
    private final String description;
    private final Consumer<String> consumer;
    private final Runnable action;

    public Command(String name, String description, Runnable action){
        this.name = name;
        this.description = description;
        this.action = action;
        this.consumer = null;
    }

    public Command(String name, String description, Consumer<String> consumer){
        this.name = name;
        this.description = description;
        this.consumer = consumer;
        this.action = null;
    }


    public String getName() {
        return name;
    }

    public String getDescription() {
        return description;
    }

    public Runnable getAction() {
        return action;
    }

    public boolean isConsumer(){
        return consumer !=null;
    }


    //For consumers
    public void execute(String args){
        if (consumer != null) {
            consumer.accept(args);
        }

    }

    //For actions
    public void execute(){
        if (action != null) {
            action.run();
        }

    }
}
