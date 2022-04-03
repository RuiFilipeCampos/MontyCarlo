import type { NextPage } from "next";

import React from "react";
import * as ch from "@chakra-ui/react";
import { useRouter } from "next/router";

const Demo: NextPage = () => {
  const [layers, setLayers] = React.useState(3);
  const handleChange = (value) => setLayers(value);

  return (
    <>
      <ch.Flex w="100vw" h="full">
        <ch.Flex shadow="md" w="full" h="full" p="10" m="50">
          <ch.VStack w="full">
            <ch.HStack w="full">
              <ch.Heading>The Onion</ch.Heading>
              <ch.Spacer />
            </ch.HStack>
            <ch.Spacer />
            <ch.FormControl>
              <ch.FormLabel>How many layers?</ch.FormLabel>
              <ch.HStack>
                <ch.NumberInput
                  step={1}
                  defaultValue={layers}
                  value={layers}
                  min={1}
                  max={10}
                  onChange={handleChange}
                >
                  <ch.NumberInputField value={layers} />
                  <ch.NumberInputStepper>
                    <ch.NumberIncrementStepper />
                    <ch.NumberDecrementStepper />
                  </ch.NumberInputStepper>
                </ch.NumberInput>

                <ch.Slider
                  flex="1"
                  focusThumbOnChange={false}
                  value={layers}
                  min={1}
                  max={10}
                  onChange={handleChange}
                >
                  <ch.SliderTrack>
                    <ch.SliderFilledTrack />
                  </ch.SliderTrack>
                  <ch.SliderThumb
                    fontSize="sm"
                    boxSize="32px"
                    children={layers}
                  />
                </ch.Slider>
              </ch.HStack>

              <ch.FormHelperText>
                Maximum of 10 layers due to performance issues. These
                performance restrictions are meant to protect the Server. In the
                full version there are no restriction.
              </ch.FormHelperText>
            </ch.FormControl>
            <FormNavigation previous="/demo" next={`/onion/${layers}/`} />
          </ch.VStack>
        </ch.Flex>
      </ch.Flex>
    </>
  );
};

export default Demo;

const LayerSelector = () => {
  const [selected, setSelected] = React.useState("none");

  const image_options = {
    none: "",
    sphere: "",
    onion:
      "https://user-images.githubusercontent.com/63464503/124515938-880a8f80-ddd8-11eb-9439-409381b5124a.png",
    "sphere-cut":
      "https://user-images.githubusercontent.com/63464503/137699819-3dad4fb6-7e76-4a5c-89b4-86924a62105c.png",
  };

  const handleSelection = (event) => setSelected(event.target.value);

  return (
    <>
      <ch.HStack w="full">
        <ch.Heading>The Onion</ch.Heading>
        <ch.Spacer />
      </ch.HStack>
      <ch.Spacer />
      <ch.FormControl>
        <ch.FormLabel>How many layers?</ch.FormLabel>
        <ch.NumberInput step={1} defaultValue={3} min={1} max={10}>
          <ch.NumberInputField />
          <ch.NumberInputStepper>
            <ch.NumberIncrementStepper />
            <ch.NumberDecrementStepper />
          </ch.NumberInputStepper>
        </ch.NumberInput>

        <ch.FormHelperText>
          Maximum of 10 layers due to performance issues. These performance
          restrictions are meant to protect the Server. In the full version
          there are no restriction.
        </ch.FormHelperText>
      </ch.FormControl>
    </>
  );
};

const FormNavigation = ({ previous, next }) => {
  const router = useRouter();

  return (
    <ch.HStack w="full" align="left">
      <ch.Spacer />
      <ch.Button onClick={() => router.push(previous)}>Previous</ch.Button>
      <ch.Button colorScheme="facebook" onClick={() => router.push(next)}>
        Next
      </ch.Button>
    </ch.HStack>
  );
};
